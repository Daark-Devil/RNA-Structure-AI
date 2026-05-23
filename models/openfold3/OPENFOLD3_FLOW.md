# OpenFold3 RNA Pipeline :- Complete Flow

How one RNA is predicted with this modified OpenFold3, and how the four
benchmark conditions are produced. Built and run on an HPC cluster (SLURM).
Paths shown are the originals - replace with your own.

---

## Environments

| Env | Purpose |
|-----|---------|
| `openfold3_clean` | Unmodified OpenFold3 (baseline reference for diffing) |
| `openfold3_rna_templ` | Modified OpenFold3 (RNA templates + RNA featurization + SS token bonds) |

The two differ by exactly **9 source files** (see `edited_files/` and `patches/`).

---

## Two-stage design

The pipeline runs in two stages. **Stage 2 reuses Stage 1's outputs** — it does
not regenerate templates or CIFs.

```
STAGE 1 : TEMPLATE PIPELINE        (run_one_rna_of3_template.sh)
  INPUT:  RNA.seq + RNA.2d
  STEPS:
    1. VFold3D motif extraction        (vfold3D_motif.o)
    2. template search                 (vfold3D_template_finder)
    3. trim + OF3 template prep        (trim_template_pdbs.py,
                                        prepare_of3_templates_from_vfold.py)
    3.5 CIF fixes                      (fix_motif_cifs_v7.py, .sto header fix,
                                        runner.yml flag fix, cache wipe)
    4. OF3 predict (templates only)    (run_openfold predict --use-templates True)
  PRODUCES (reused later):
    of3_structures/*.cif               <- the template CIFs (auto-generated here)
    RNA_templates.sto                  <- template alignment
    RNA_query_with_templates.json      <- template chain ids

           |  Stage 2 copies of3_structures/*.cif + RNA_templates.sto
           v
STAGE 2 : MSA / TOKEN-BOND LAYER   (run_4cond_batch.sh)
  REUSES:  of3_structures/*.cif, RNA_templates.sto  (copied, identical)
  GENERATES: synthetic MSA            (ss_to_synthetic_msa.py)
  BUILDS:  per-condition JSON + runner.yml
  RUNS OF3 in 4 conditions (below)
```

The CIF being identical across stages is confirmed: the Stage-2 CIF is a byte-for-
byte copy of the Stage-1 CIF. **No manual CIF creation at any point** - your prep
scripts make it in Stage 1, Stage 2 reuses it.

---

## The 4 benchmark conditions (Stage 2)

| # | Name | use_msas | SS token bonds env | Meaning |
|---|------|----------|--------------------|---------|
| 1 | Baseline | false | no | templates only |
| 2 | MSA | true (+ main_msa_file_paths) | no | templates + synthetic MSA |
| 3 | Token bonds | false | `SS_TOKEN_BONDS_FILE` set | templates + SS token bonds |
| 4 | Combined | true (+ MSA path) | `SS_TOKEN_BONDS_FILE` set | templates + MSA + token bonds |

Two switches control everything:
- **Synthetic MSA**: in the query JSON — `use_msas: true` + `main_msa_file_paths`
  pointing at the synthetic MSA folder. No code edit needed.
- **SS token bonds**: environment variables read by the patched featurization code:
  ```bash
  export SS_TOKEN_BONDS_FILE=<RNA>.2d
  export SS_TOKEN_BONDS_CHAIN_LEN=<length>
  ```

---

## How inputs reach the model (no raw CLI files)

`run_openfold predict` takes only a **query JSON** and a **runner YAML**. Everything
else is referenced by path inside those:

```
run_openfold predict \
  --query-json   RNA_with_synthetic_msa.json   # contains: sequence, MSA path,
                                                #           .sto path, template ids
  --runner-yaml  runner_RNA.yml                 # contains: structure_directory
                                                #           (where CIFs live) + cache
  --inference-ckpt-path  of3-p2-145k.pt
  --use-msa-server False
  --use-templates True
  --num-diffusion-samples 10
  --num-model-seeds 3
  --output-dir   ./of3_output_msa
```

- **MSA** -> JSON `main_msa_file_paths` (a folder of .sto)
- **Template alignment** -> JSON `template_alignment_file_path` (the .sto)
- **Template CIFs** -> runner YAML `structure_directory`; OF3 builds its own
  precache (`create_precache: true`) into `template_data/template_cache/`
- **SS token bonds** -> environment variables (read by patched code)

So CIFs are NOT a CLI argument — they are found via the runner YAML's
`structure_directory`, and matched to the `.sto` chain ids (e.g. `m0_A`).

---

## Synthetic MSA generator

`ss_to_synthetic_msa.py` (standalone, imports only argparse/random/pathlib — no
OpenFold3 dependency). Turns a dot-bracket SS into a Stockholm MSA using
RNA pairing biology:
- pair preferences (A-U 85%, G-C 85%, G-U wobble)
- length-tier sequence count (256–768 by length)
- stem-length / stem-position / loop-type aware mutation rates

```bash
python ss_to_synthetic_msa.py \
    --seq RNA.seq --ss RNA.2d \
    --out synthetic_msa/chain_A/rfam_hits.sto \
    --query_name RNA_A
```

---

## Outputs

Per seed, per diffusion sample: `*_model.cif`. RMSD vs native (C1' Kabsch,
`compare_rna_kabsch.py`) selects best-of-N. Per-RNA and aggregate summaries:
`all_rnas_rmsd_sorted.tsv` (PIPELINE_A), `all_rnas_rmsd_sorted_PIPELINE_B.tsv`.

---

## Reproduce one RNA

```bash
conda activate openfold3_rna_templ

# Stage 1: template pipeline (produces CIF + .sto)
bash run_one_rna_of3_template.sh 1Z30

# Stage 2: add synthetic MSA + run (reuses Stage 1 CIF/.sto)
bash run_ss_msa_pipeline.sh 1Z30 \
     <path>/1Z30.seq <path>/1Z30.2d \
     <path>/batch_runs/withtemplate/1Z30 \
     <path>/1Z30.pdb     # optional, for RMSD
```

Paths to change for another machine: project root, rna_workdirs, checkpoint
(`params/of3-p2-145k.pt`), and the VFold3D binary locations.

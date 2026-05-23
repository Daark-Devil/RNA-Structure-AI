# OpenFold3 :- RNA Template + Synthetic-MSA Pipeline

Predicting RNA 3D structure with OpenFold3, modified to accept RNA templates and
(optionally) a synthetic secondary-structure-derived MSA.

The pipeline runs in two stages, two commands:

```
Stage 1  bash run_one_rna_of3_template.sh <RNA>          # templates -> makes CIF + .sto
Stage 2  bash run_ss_msa_pipeline.sh <RNA> <seq> <ss> <stage1_dir> [pdb]   # + synthetic MSA
```

Stage 2 reuses the CIF and .sto produced by Stage 1 (it does not regenerate them).


## Folder map (what is where)

```
edited_files/   the 9 modified OpenFold3 source files (+ README on how to replace)
patches/        the same 9 changes as unified diffs
scripts/        run_one_rna_of3_template.sh   (Stage 1)
                run_ss_msa_pipeline.sh        (Stage 2)
                prepare_of3_templates_from_vfold.py
                fix_motif_cifs_v7.py
                trim_template_pdbs.py
                ss_to_synthetic_msa.py        (synthetic MSA generator)
                compare_rna_kabsch.py         (RMSD)
                templates/  with_msa.json.tpl, template_only.json.tpl, runner.yml.tpl
example_17RA/   one complete run: inputs -> intermediates -> output
results/        all_rnas_rmsd_sorted.tsv (+ PIPELINE_B baseline)
docs/           OPENFOLD3_FLOW.md (full flow explanation)
```


## Prerequisites
- Official OpenFold3 installed, then the 9 files replaced (see `edited_files/README.md`)
- OpenFold3 checkpoint `of3-p2-145k.pt`
- VFold3D binaries (`vfold3D_motif.o`, `vfold3D_template_finder`) + `helix.pdb`
- A conda env (here `openfold3_rna_templ`) and CUDA toolkit module


## STAGE 1 : Template pipeline

```bash
conda activate openfold3_rna_templ
bash run_one_rna_of3_template.sh <RNA>
```

Argument:
- `<RNA>` : RNA name, e.g. `1Z30`. The script looks for `<RNA>.seq` and `<RNA>.2d`
  inside the workdir base (see paths to change below).

Produces (under `<BATCH_BASE>/withtemplate/<RNA>/`):
- `of3_structures/*.cif`            template CIFs (auto-generated here)
- `<RNA>_templates.sto`             template alignment
- `<RNA>_query_with_templates.json` template chain ids
- `output_of3-p2-145k/`             template-only prediction


## STAGE 2 : Add synthetic MSA and run

```bash
bash run_ss_msa_pipeline.sh <RNA> <SEQ_FILE> <SS_FILE> <STAGE1_DIR> [REF_PDB]
```

Arguments (in order):
- `<RNA>`        RNA name, e.g. `1Z30`
- `<SEQ_FILE>`   path to `<RNA>.seq` (one-line sequence)
- `<SS_FILE>`    path to `<RNA>.2d` (two lines: sequence, then dot-bracket)
- `<STAGE1_DIR>` the Stage-1 output dir for this RNA
                 (e.g. `.../batch_runs/withtemplate/1Z30`) — this is where the
                 script copies the CIF + .sto from
- `[REF_PDB]`    optional native `.pdb`; if given, prints best/mean RMSD

Produces (under `~/ss_msa_test/<RNA>/`):
- `synthetic_msa/chain_A/rfam_hits.sto`  the generated synthetic MSA
- `<RNA>_with_synthetic_msa.json`        query with `use_msas: true`
- `<RNA>_baseline.json`                  query with `use_msas: false`
- `of3_output_msa/.../*_model.cif`       predicted structures


## PATHS YOU MUST CHANGE

### In `run_one_rna_of3_template.sh`
| Line | Variable | Change to |
|------|----------|-----------|
| 54 | `WORKDIR_BASE=` | your folder of `<RNA>/<RNA>.seq` + `.2d` inputs |
| 55 | `BATCH_BASE=`    | where Stage-1 outputs should be written |
| 57 | `VFROOT=`        | your VFold3D `Data` directory |
| 60 | `TRIM_PY=`       | path to `trim_template_pdbs.py` |
| 62 | `PREP_PY=`       | path to `prepare_of3_templates_from_vfold.py` |
| 64 | `CKPT=`          | path to your `of3-p2-145k.pt` checkpoint |
| 200 | hardcoded path to `fix_motif_cifs_v7.py` | your path to that script |
| 30-32 | `module load cuda...` | your cluster's CUDA module name |

### In `run_ss_msa_pipeline.sh`
| Line | Variable | Change to |
|------|----------|-----------|
| 26 | `PIPELINE=` | path to the `ss_msa_pipeline` folder (holds `templates/`) |
| 27 | `WORKDIR=`  | where Stage-2 outputs should be written |
| 28 | `GEN=`      | path to `ss_to_synthetic_msa.py` |
| 29 | `CKPT=`     | path to `of3-p2-145k.pt` |
| 104 | `KABSCH=`  | path to `compare_rna_kabsch.py` |
| 33 | `module load cuda...` | your cluster's CUDA module name |

Note: Stage 2 reads its JSON/YAML templates from `$PIPELINE/templates/`
(`with_msa.json.tpl`, `template_only.json.tpl`, `runner.yml.tpl`), so make sure
`PIPELINE` points at the folder that contains `templates/`.


## Synthetic MSA generator (standalone)

`ss_to_synthetic_msa.py` turns SS into a Stockholm MSA using RNA pairing biology
(A-U/G-C ~85%, G-U wobble; stem/loop-aware mutation rates; sequence count scales
with length). No OpenFold3 dependency. Stage 2 calls it for you, or run directly:

```bash
python ss_to_synthetic_msa.py \
    --seq <RNA>.seq --ss <RNA>.2d \
    --out synthetic_msa/chain_A/rfam_hits.sto \
    --query_name <RNA>_A
```


## How inputs reach the model

`run_openfold predict` takes only the query JSON and runner YAML. Everything else
is referenced by path inside them:
- synthetic MSA  -> JSON `main_msa_file_paths`
- template .sto  -> JSON `template_alignment_file_path`
- template CIFs  -> runner YAML `structure_directory` (OF3 builds its own precache)

CIFs are never passed as a CLI argument.

# RNA-Structure-AI

Improving RNA 3D structure prediction by adding **template** and **secondary-
structure** guidance to open-source AlphaFold3 reproductions. Three model families
are explored; each is modified to accept RNA structural guidance that the stock
code does not support, then benchmarked against the unmodified baseline.

## Models

| Model | Approach | Status |
|-------|----------|--------|
| **OpenFold3** | RNA templates + synthetic SS-MSA (9 edited files) | included |
| **Boltz-2** | RNA templates + SS contact constraints (2 edited files) | included |
| **GraphaRNA** | GNN + diffusion with N motif templates | planned |

All work is built on the open-source projects and shared as edits/patches on top
of the official code :- install the upstream model, then apply the changes here.

## Repository layout

```
models/
  openfold3/   RNA template + synthetic-MSA pipeline
    edited_files/  patches/  scripts/  example_17RA/  example_1Z30_msa/
    results/  docs/  README.md
  boltz2/      SS constraints + RNA templates
    edited_files/  patches/  scripts/  example_ss_constraints/
    example_template/  results/  README.md
figures/       flow + structure diagrams, result visualizations
```

Each model folder is self-contained: edited source files (with a README on how to
apply them to a fresh install), the pipeline scripts, a complete worked example
(input -> intermediates -> output), and benchmark results.

## What was changed and why

Stock AlphaFold3 reproductions support protein templates and protein MSAs, but not
RNA equivalents. This work adds:

- **OpenFold3** : accepts RNA 3D templates (C1' distances, P-C4'-C1' frames) and an
  optional synthetic MSA generated from the RNA secondary structure, giving the
  model coevolution-like signal without a real alignment.
- **Boltz-2** : accepts RNA template chains (frames from C4'-C1'-glycosidic N), and
  a generator that turns dot-bracket SS into contact constraints (works on stock
  Boltz, no patch needed).

Structures are evaluated by C1' Kabsch RMSD vs the native, best-of-N over diffusion
samples and seeds.

## Quick start (one command per model)

Each model's full pipeline runs from a single entry command. See the model README
for arguments and the paths to change.

**OpenFold3** (template stage, then synthetic-MSA stage):
```bash
conda activate openfold3_rna_templ
bash models/openfold3/scripts/run_one_rna_of3_template.sh 1Z30
bash models/openfold3/scripts/run_ss_msa_pipeline.sh \
     1Z30 1Z30.seq 1Z30.2d <stage1_dir> 1Z30.pdb
```

**Boltz-2** (pick an approach):
```bash
conda activate boltz2
# SS constraints (stock Boltz-2)
python models/boltz2/scripts/generate_boltz_yaml.py \
       --rna 1A1T --seq 1A1T.seq --ss 1A1T.2d --out 1A1T.yaml
boltz predict 1A1T.yaml --diffusion_samples 5 --out_dir results_1A1T --override

# template (patched Boltz-2)
python models/boltz2/scripts/generate_template_yaml.py \
       --rna 1A1T --pdb 1A1T.pdb --out 1A1T_template.yaml --cif_out 1A1T.cif
boltz predict 1A1T_template.yaml --out_dir results_1A1T_template --override
```

## Results summary

Across a 22-RNA benchmark, OpenFold3 with synthetic MSA was the strongest single
method overall, with large gains on several RNAs (e.g. 1P5P 22.45 -> 3.54 A,
1E7K 8.98 -> 0.94 A). Boltz-2 templates won on some RNAs (e.g. 3G8T 2.80 A),
showing the two approaches are complementary rather than one strictly dominating.
Per-RNA tables and head-to-head comparisons are in each model's `results/` folder.
Some RNAs did not improve (e.g. 1U9S with synthetic MSA) :- these are reported too.

## Notes
- No proprietary weights or large datasets are included; install upstream models
  for those.
- Paths in scripts use placeholders / are passed as arguments :- see each model
  README for the specific paths to set.

# Boltz-2 : Step-by-Step Workflow (start to finish)

Complete flow from applying the code edits to getting a prediction, scoring it,
and visualizing it. Two approaches are shown; pick one per run.

---

## Step 0 : Install + (for templates) apply the edits

Install official Boltz-2 into a conda env. For the **template** approach only,
apply the 2 edited files (the SS-constraint approach needs no edits):

```bash
BOLTZ=$(python -c "import boltz, os; print(os.path.dirname(boltz.__file__))")

# back up, then replace
cp "$BOLTZ/data/parse/schema.py"     "$BOLTZ/data/parse/schema.py.orig"
cp "$BOLTZ/data/tokenize/boltz2.py"  "$BOLTZ/data/tokenize/boltz2.py.orig"
cp -r edited_files/* "$BOLTZ/"

python -c "import boltz; print('import ok')"
```
Full detail + patch alternative: `edited_files/README.md`.

---

## Approach 1 : SS constraints (stock Boltz-2)

### Step 1 : inputs
A one-line `.seq` and a `.2d` (sequence + dot-bracket), e.g. `1A1T.seq`, `1A1T.2d`.

### Step 2 : generate the YAML
```bash
python scripts/generate_boltz_yaml.py \
    --rna 1A1T --seq 1A1T.seq --ss 1A1T.2d --out 1A1T.yaml
```
Turns each base pair into a Boltz contact constraint. (See `scripts/README.md`
for the YAML shape.)

### Step 3 : predict
```bash
boltz predict 1A1T.yaml \
    --diffusion_samples 5 --recycling_steps 3 --sampling_steps 200 \
    --out_dir results_1A1T --override
```

### Step 4 : output
```
results_1A1T/boltz_results_1A1T/predictions/1A1T/1A1T_model_0.cif   (+ model_1..4)
```

---

## Approach 2 : template (patched Boltz-2)

### Step 1 : input
An RNA template structure as PDB (e.g. `1A1T.pdb`).

### Step 2 : generate the YAML (also writes the CIF it points to)
```bash
python scripts/generate_template_yaml.py \
    --rna 1A1T --pdb 1A1T.pdb \
    --out 1A1T_template.yaml --cif_out 1A1T.cif
```
Note: gemmi renames chain A to `Axp` during PDB->CIF; the generator sets
`template_id: [Axp]` automatically — do not change it.

### Step 3 : predict (best-of-5 seeds)
```bash
for S in 0 1 2 3 4; do
  boltz predict 1A1T_template.yaml --seed $S --out_dir results_1A1T_s$S --override
done
```

### Step 4 : output
```
results_1A1T_s<seed>/boltz_results_*/predictions/*/*_model_0.cif
```

---

## Step 5 : score (RMSD vs native)

Use the C1' Kabsch tool (from the OpenFold3 section) to compare a prediction to
the native structure:
```bash
python compare_rna_kabsch.py --ref 1A1T.pdb --pred <predicted>_model_0.cif
# prints C1'_KABSCH_RMSD=<value>
```
Take best-of-N over the samples/seeds.

---

## Step 6 : visualize

3-panel image (prediction / native / overlay), aligned on C1':
```bash
python scripts/viz_boltz_comparison.py \
    --rna 3G8T \
    --pred   <predicted>_model_0.cif \
    --native 3G8T.pdb \
    --out-dir viz_out \
    --rmsd-tool compare_rna_kabsch.py
```
Produces `viz_out/3G8T/boltz_vs_native.png`. Requires PyMOL. Example output is in
`results/figures/`.

---

## Summary of the flow

```
(template only) apply 2 edits
        |
inputs (.seq/.2d  or  .pdb)
        |
generate YAML  (scripts/generate_*.py)
        |
boltz predict  -> *_model_0.cif
        |
compare_rna_kabsch.py  -> RMSD (best-of-N)
        |
viz_boltz_comparison.py -> 3-panel PNG
```

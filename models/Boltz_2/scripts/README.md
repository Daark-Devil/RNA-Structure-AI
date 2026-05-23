# Boltz-2 YAML Generators

Two generators produce the two kinds of Boltz-2 input YAML used here. Both take a
single RNA and write a ready-to-run YAML.

| Generator | Approach | Input | Needs patched Boltz? |
|-----------|----------|-------|----------------------|
| `generate_boltz_yaml.py` | SS constraints | `.seq` + `.2d` | No |
| `generate_template_yaml.py` | Template | `.pdb` | Yes |

---

## 1. generate_boltz_yaml.py  (SS constraints)

Parses a dot-bracket secondary structure into base pairs and writes each pair as a
Boltz **contact constraint** between the two paired residues.

### Inputs
- `--seq` one-line sequence file, e.g. `GGACUAGCGGAGGCUAGUCC`
- `--ss`  dot-bracket file (two lines: sequence, then structure), e.g.
  ```
  GGACUAGCGGAGGCUAGUCC
  ((((((((....))))))))
  ```

### Run
```bash
python generate_boltz_yaml.py \
    --rna 1A1T --seq 1A1T.seq --ss 1A1T.2d --out 1A1T_ss.yaml
```

### Output YAML (shape)
```yaml
version: 1
sequences:
  - rna:
      id: A
      sequence: GGACUAGCGGAGGCUAGUCC
constraints:
  - contact:
      token1: [A, 1]
      token2: [A, 20]
      max_distance: 4.0
      force: true
  - contact:
      token1: [A, 2]
      token2: [A, 19]
      ...
```
Each base pair (i, j) from the dot-bracket becomes one `contact` constraint with
`max_distance: 4.0` and `force: true`. Supports nested/pseudoknot brackets
( ), [ ], { }, < >.

### Predict
```bash
boltz predict 1A1T_ss.yaml \
    --diffusion_samples 5 --recycling_steps 3 --sampling_steps 200 \
    --out_dir results_1A1T --override
```

---

## 2. generate_template_yaml.py  (template)

Converts an RNA template PDB to CIF (via gemmi) and writes a YAML referencing it.
Requires the patched Boltz-2.

### Inputs
- `--pdb` the template structure (PDB)
- `--out` output YAML path
- `--cif_out` where to write the converted CIF (the YAML points at this)

### Run
```bash
python generate_template_yaml.py \
    --rna 1A1T --pdb 1A1T.pdb \
    --out 1A1T_template.yaml --cif_out 1A1T.cif
```

### Output YAML (shape)
```yaml
version: 1
sequences:
  - rna:
      id: A
      sequence: GGACUAGCGGAGGCUAGUCC
templates:
  - cif: ./1A1T.cif
    chain_id: [A]
    template_id: [Axp]
```

### The Axp detail
gemmi renames chain A to **Axp** during PDB->CIF conversion, so `template_id` must
be `[Axp]`. The generator sets this automatically from the converted CIF — do not
change it back to A.

### Predict (best-of-5 seeds)
```bash
for S in 0 1 2 3 4; do
  boltz predict 1A1T_template.yaml --seed $S --out_dir results_1A1T_s$S --override
done
```

---

## Output location (both approaches)
```
results_*/boltz_results_*/predictions/<name>/<name>_model_0.cif   (+ _model_1..N)
```
Pick best-of-N by RMSD vs a reference (e.g. with the C1' Kabsch tool from the
OpenFold3 section).

## Notes
- No MSA is used in either approach.
- SS constraints work on stock Boltz-2; templates require the 2 edited files
  (see `../edited_files/README.md`).

---

## 3. viz_boltz_comparison.py  (visualization)

3-panel comparison PNG: Boltz prediction / native / overlay (aligned on C1').
Requires PyMOL.

```bash
python viz_boltz_comparison.py \
    --rna 3G8T \
    --pred   <predicted>_model_0.cif \
    --native 3G8T.pdb \
    --out-dir viz_out \
    --rmsd-tool /path/to/compare_rna_kabsch.py
```
Output: `viz_out/<RNA>/boltz_vs_native.png`. Overlay colors: native black,
prediction green.

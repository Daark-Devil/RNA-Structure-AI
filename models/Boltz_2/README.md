# Boltz-2 : RNA Structure Prediction (SS Constraints + Templates)

Two ways to guide Boltz-2 RNA structure prediction:

1. **SS constraints** :- convert a secondary structure (dot-bracket) into base-pair
   contact constraints in the Boltz YAML. Works on **stock Boltz-2**.
2. **Templates** :- supply an RNA 3D template (CIF) and have Boltz-2 use its
   geometry. Requires a **patched Boltz-2** (2 edited files; stock Boltz rejects
   RNA template chains).

No MSA is used in either approach.

## Folder map

```
README.md                    this file
HOW_TO_RUN.txt               quick command reference
edited_files/                the 2 modified Boltz source files + how to apply
  README.md                  find package path, replace or patch, verify
  data/parse/schema.py
  data/tokenize/boltz2.py
patches/                     the 2 changes as unified diffs
scripts/
  generate_boltz_yaml.py     SS-constraint YAML generator (approach 1)
  generate_template_yaml.py  template YAML generator (approach 2)
example_ss_constraints/      1A1T: seq + 2d -> yaml -> predicted cif
example_template/            1A1T: template cif -> yaml -> predicted cif
results/                     benchmark tables
```

## The code changes (approach 2 only)

| File | Lines | Change |
|------|-------|--------|
| `data/parse/schema.py` | 16 | Accept RNA/DNA template chains (was protein-only) |
| `data/tokenize/boltz2.py` | 33 | Build RNA/DNA template frames from C4'->C1'->glycosidic N (N9 purines, N1 pyrimidines), analogous to the protein N-CA-C frame |

See `edited_files/README.md` to apply them. Approach 1 (SS constraints) needs no
code changes.

## Approach 1 : SS constraints (stock Boltz-2)

Inputs: a one-line `.seq` and a `.2d` (sequence + dot-bracket).

```bash
conda activate boltz2
python scripts/generate_boltz_yaml.py \
    --rna 1A1T --seq 1A1T.seq --ss 1A1T.2d --out 1A1T.yaml

boltz predict 1A1T.yaml \
    --diffusion_samples 5 --recycling_steps 3 --sampling_steps 200 \
    --out_dir results_1A1T --override
```
Output: `results_1A1T/boltz_results_1A1T/predictions/1A1T/1A1T_model_0.cif` (+ 1..4).

The generator parses the dot-bracket into base pairs and writes them as Boltz
contact constraints. Example in `example_ss_constraints/`.

## Approach 2 : Template (patched Boltz-2)

Input: an RNA template structure as PDB. The generator converts it to CIF and
writes a YAML that references it.

```bash
conda activate boltz2     # with the patched Boltz-2 installed
python scripts/generate_template_yaml.py \
    --rna 1A1T --pdb 1A1T.pdb \
    --out 1A1T_template.yaml --cif_out 1A1T.cif

# 5 seeds for best-of-5
for S in 0 1 2 3 4; do
  boltz predict 1A1T_template.yaml --seed $S --out_dir results_1A1T_s$S --override
done
```
Output: `results_1A1T_s<seed>/boltz_results_*/predictions/*/*_model_0.cif`.

**Proof of concept note:** stock Boltz-2 rejects RNA template chains. After
patching, Boltz-2 accepts an RNA template CIF and builds frames for it. This was
validated on 1A1T by supplying a reference RNA structure as the template and
confirming the run completes and uses the template. In production, the template
would be a found/homologous structure (e.g. from the VFold motif pipeline used in
the OpenFold3 section) rather than the reference itself. Example in
`example_template/` (the included `1A1T.cif` is the reference-as-template).

## The Axp gotcha (template approach)

When `gemmi` converts PDB->CIF, chain A is renamed to **Axp**. The template YAML
must therefore use `template_id: [Axp]`. The generator sets this automatically —
do not change it back to A.

## Paths to change
The generators take all paths as arguments, so there is nothing hardcoded to edit.
For the template example YAML, the `cif:` field points to `./1A1T.cif` — point it
at your own template CIF (or the one the generator writes via `--cif_out`).

## Environment
```bash
conda activate boltz2
# approach 2 also requires the 2 edited files applied (see edited_files/README.md)
```

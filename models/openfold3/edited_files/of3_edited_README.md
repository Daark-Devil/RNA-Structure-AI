# OpenFold3 — RNA Template Support (Modified Files)

These are the **9 OpenFold3 source files** modified to add RNA template support
(stock OpenFold3 supports protein templates only). This folder mirrors the internal
OpenFold3 package structure, so replacing the originals is a single recursive copy.

> Base project: OpenFold3 (open-source AlphaFold3 reproduction). Install the official
> OpenFold3 first, then replace these 9 files. Modifications by Devansh Pancholi.

## The 9 modified files

| File (path inside `openfold3/`) | What changed |
|---|---|
| `core/data/pipelines/sample_processing/template.py` | Accept RNA chains in template chain filter (was protein-only) |
| `core/data/io/sequence/template.py` | Accept `mol:rna` headers in template alignment parsing |
| `core/data/primitives/featurization/template.py` | RNA (A,G,C,U): use C1' for distances, P-C4'-C1' for frames |
| `core/data/primitives/featurization/structure.py` | RNA-aware structure / token-bond featurization |
| `core/data/primitives/structure/template.py` | RNA template structure handling |
| `core/data/pipelines/preprocessing/template.py` | RNA template atom QC (P, C4', C1' instead of N, CA, C, CB) |
| `core/data/primitives/quality_control/asserts.py` | Relax protein-only asserts so RNA passes QC |
| `core/data/framework/single_datasets/inference.py` | RNA handling in inference dataset |
| `projects/of3_all_atom/config/model_config.py` | Attention-kernel toggle (DeepSpeed off / LMA on); does NOT change MSA/template logic |

## How to install (replace the originals)

The internal structure under `site-packages/openfold3/` is identical on every
machine — only the prefix path differs. So we find the prefix once, then copy.

### Step 1 — find your OpenFold3 install path
```bash
OF3=$(python -c "import openfold3, os; print(os.path.dirname(openfold3.__file__))")
echo "$OF3"
# e.g. /opt/conda/envs/yourenv/lib/python3.11/site-packages/openfold3
```

### Step 2 — back up the 9 originals (so you can revert)
```bash
for f in \
  core/data/framework/single_datasets/inference.py \
  core/data/io/sequence/template.py \
  core/data/pipelines/preprocessing/template.py \
  core/data/pipelines/sample_processing/template.py \
  core/data/primitives/featurization/structure.py \
  core/data/primitives/featurization/template.py \
  core/data/primitives/quality_control/asserts.py \
  core/data/primitives/structure/template.py \
  projects/of3_all_atom/config/model_config.py
do
  cp "$OF3/$f" "$OF3/$f.orig"
done
echo "backups created (.orig)"
```

### Step 3 — replace with the modified files
This folder mirrors the OpenFold3 tree, so one recursive copy puts every file in
the right place:
```bash
cp -r edited_files/* "$OF3/"
```

### Step 4 — verify it still imports
```bash
python -c "import openfold3; print('import ok')"
```

## How to revert
```bash
# restore each original from its .orig backup
for f in \
  core/data/framework/single_datasets/inference.py \
  core/data/io/sequence/template.py \
  core/data/pipelines/preprocessing/template.py \
  core/data/pipelines/sample_processing/template.py \
  core/data/primitives/featurization/structure.py \
  core/data/primitives/featurization/template.py \
  core/data/primitives/quality_control/asserts.py \
  core/data/primitives/structure/template.py \
  projects/of3_all_atom/config/model_config.py
do
  cp "$OF3/$f.orig" "$OF3/$f"
done
```

## Notes
- The unified diffs of these same changes are in `../patches/` if you only want to
  see what changed rather than replace whole files.
- Tested with OpenFold3 on Python 3.11.

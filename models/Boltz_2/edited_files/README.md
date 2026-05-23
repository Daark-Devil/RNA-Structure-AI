# Boltz-2 : RNA Template Support (Modified Files)

These are the **2 Boltz-2 source files** modified to add RNA (and DNA) template
support. Stock Boltz-2 only accepts protein templates; these changes let it accept
RNA template chains and build geometric frames for them.

> Base project: Boltz-2 (open-source AlphaFold3 reproduction). Install the official
> Boltz-2 first, then apply these changes. Only the template approach needs them —
> the SS-constraint approach works on stock Boltz-2.

## The 2 modified files

| File (path inside `boltz/`) | Lines | What changed |
|---|---|---|
| `data/parse/schema.py` | 16 | Accept RNA/DNA template chains (was protein-only) in two filters: the templateable-chain set and the template-chain detection |
| `data/tokenize/boltz2.py` | 33 | Build a coordinate frame for RNA/DNA template residues using **C4' -> C1' -> glycosidic N** (N9 for purines A/G, N1 for pyrimidines C/U), analogous to the protein N-CA-C frame |

Net effect: Boltz-2 accepts an RNA 3D template (as a CIF) and uses its geometry,
which the stock code only does for proteins.

## How to apply :- two ways (pick one)

First, find your Boltz-2 install path:
```bash
BOLTZ=$(python -c "import boltz, os; print(os.path.dirname(boltz.__file__))")
echo "$BOLTZ"
# e.g. /opt/conda/envs/yourenv/lib/python3.10/site-packages/boltz
```

### Way A :- replace whole files (simplest)
This folder mirrors the Boltz tree, so one recursive copy puts both files in place.
Back up the originals first:
```bash
cp "$BOLTZ/data/parse/schema.py"      "$BOLTZ/data/parse/schema.py.orig"
cp "$BOLTZ/data/tokenize/boltz2.py"   "$BOLTZ/data/tokenize/boltz2.py.orig"

cp -r edited_files/* "$BOLTZ/"
```

### Way B :- apply patches (only changes the edited lines)
```bash
cp "$BOLTZ/data/parse/schema.py"      "$BOLTZ/data/parse/schema.py.orig"
cp "$BOLTZ/data/tokenize/boltz2.py"   "$BOLTZ/data/tokenize/boltz2.py.orig"

patch "$BOLTZ/data/parse/schema.py"     < patches/schema.py.patch
patch "$BOLTZ/data/tokenize/boltz2.py"  < patches/boltz2_tokenize.py.patch
```

### Verify (either way)
```bash
python -c "import boltz; print('import ok')"
```

## Revert
```bash
cp "$BOLTZ/data/parse/schema.py.orig"     "$BOLTZ/data/parse/schema.py"
cp "$BOLTZ/data/tokenize/boltz2.py.orig"  "$BOLTZ/data/tokenize/boltz2.py"
```

## Notes
- The two `.patch` files in `../patches/` are the same changes as unified diffs.
- There are TWO files named `boltz2.py` in Boltz; the one edited is the
  **tokenizer** at `data/tokenize/boltz2.py`, NOT `model/models/boltz2.py`.
- Tested with Boltz-2 on Python 3.10.

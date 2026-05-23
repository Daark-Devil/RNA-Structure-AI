#!/usr/bin/env bash
# =============================================================================
# Script 1 — per-RNA OF3 prediction, template-only, 145k checkpoint
# Usage: bash run_one_rna_of3_template.sh <RNA>
#
# Exit codes:
#   0  = OF3 prediction completed
#   1  = bad usage / missing inputs
#   2  = length filter (RNA too long)
#   3  = template prep failed (no usable VFold templates)  -> master should skip post-processing
#   4  = OF3 prediction failed
# =============================================================================
set -u

RNA="${1:-}"
if [[ -z "$RNA" ]]; then
  echo "[ERR] usage: $0 <RNA_NAME>"
  exit 1
fi

# ----------- env setup -----------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate openfold3_rna_templ

# Load CUDA toolkit module (DeepSpeed needs nvcc at import time)
# Clear any stale CUDA_HOME first
unset CUDA_HOME

if command -v module >/dev/null 2>&1; then
  module load cuda12.8/toolkit/12.8.1 2>/dev/null || \
  module load cuda/12.8 2>/dev/null || \
  module load cuda 2>/dev/null || true
fi

if command -v nvcc >/dev/null 2>&1; then
  export CUDA_HOME="$(dirname "$(dirname "$(command -v nvcc)")")"
  export LD_LIBRARY_PATH="$CUDA_HOME/lib64:$CUDA_HOME/targets/x86_64-linux/lib:${LD_LIBRARY_PATH:-}"
else
  echo "[ERR] nvcc not found after module load — DeepSpeed will fail. Run 'module load cuda12.8/toolkit/12.8.1' manually first."
  exit 1
fi

export CUTLASS_PATH="$(python - << 'PY'
import pathlib
try:
    import cutlass_library
    print(pathlib.Path(cutlass_library.__file__).resolve().parent.joinpath("source"))
except Exception:
    print("")
PY
)"

# ----------- paths -----------
WORKDIR_BASE="/home/dpancholi1/Projects/GraphRNA_Modified/training_data/rna_workdirs"
BATCH_BASE="/home/dpancholi1/Projects/Openfold_template/batch_runs"

VFROOT="/home/dpancholi1/Projects/hou_new/Extention_Database/VfoldPipeline_standalone/Vfold3DLA/Data"
MOTIF_BIN="$VFROOT/vfold3D/vfold3D_motif.o"
TEMPLATE_BIN="$VFROOT/vfold3D/vfold3D_template_finder"
TRIM_PY="/home/dpancholi1/Projects/hou_new/Extention_Database/SAFE_scripts/trim_template_pdbs.py"
HELIX_PDB="$VFROOT/data/vfold3D/helix.pdb"
PREP_PY="/home/dpancholi1/Projects/Openfold_template/prepare_of3_templates_from_vfold.py"

CKPT="/home/dpancholi1/Projects/Openfold_template/params/of3-p2-145k.pt"
CKBASE="$(basename "$CKPT" .pt)"

NUM_DIFF=10
NUM_SEEDS=5

# ----------- inputs -----------
RNADIR="$WORKDIR_BASE/$RNA"
SEQFILE="$RNADIR/$RNA.seq"
SSFILE="$RNADIR/$RNA.2d"

if [[ ! -f "$SEQFILE" || ! -f "$SSFILE" ]]; then
  echo "[ERR] $RNA: missing seq or 2d file"
  exit 1
fi

# Length filter (<200 nt)
LEN=$(tr -d '\n' < "$SEQFILE" | wc -c)
if [ "$LEN" -ge 200 ]; then
  echo "[SKIP] $RNA length=$LEN (>=200)"
  exit 2
fi

echo "=============================="
echo "[RNA] $RNA  length=$LEN"
echo "=============================="

# ----------- output dirs -----------
TMP_BASE="$BATCH_BASE/withtemplate/$RNA"
mkdir -p "$TMP_BASE"

STATUS="$BATCH_BASE/status.tsv"
mkdir -p "$BATCH_BASE"
[[ -f "$STATUS" ]] || echo -e "rna\tstage\tstatus\tnote" > "$STATUS"

SEQ="$(head -1 "$SEQFILE")"
SS="$(tail -1 "$SSFILE")"

# ===== Step 1: motif extraction =====
printf "%s\n%s\n" "$SEQ" "$SS" > "$TMP_BASE/${RNA}_in2d.txt"
rm -f "$TMP_BASE/${RNA}_vfold3D.err" "$TMP_BASE/${RNA}_template_search.err"
rm -rf "$TMP_BASE/templates" "$TMP_BASE/of3_templates" "$TMP_BASE/of3_templates_cif" "$TMP_BASE/of3_structures" "$TMP_BASE/template_data"
touch "$TMP_BASE/excluded_pdblist.txt"

"$MOTIF_BIN" "$TMP_BASE/${RNA}_in2d.txt" "$TMP_BASE/${RNA}_vfold3D.mtf" "$TMP_BASE/${RNA}_vfold3D.err" 1

if [[ ! -s "$TMP_BASE/${RNA}_vfold3D.mtf" ]]; then
  echo -e "$RNA\tmotif_extract\tfail\tno mtf" >> "$STATUS"
  echo "[FAIL] $RNA: motif extraction produced no .mtf"
  exit 3
fi

awk 'BEGIN{keep=1} /^#####/{keep=0} keep{print}' "$TMP_BASE/${RNA}_vfold3D.mtf" > "$TMP_BASE/${RNA}_vfold3D_clean.mtf"
echo -e "$RNA\tmotif_extract\tok\tmtf created" >> "$STATUS"

# ===== Step 2: template search =====
(
  cd "$TMP_BASE" || exit 1
  "$TEMPLATE_BIN" "$VFROOT" "$SEQFILE" "${RNA}_vfold3D_clean.mtf" "${RNA}_template_search.err" 1 excluded_pdblist.txt
)

if [[ ! -f "$TMP_BASE/${RNA}_vfold3D_clean.mtf.templates.tsv" ]]; then
  echo -e "$RNA\ttemplate_search\tfail\tno templates.tsv" >> "$STATUS"
  echo "[FAIL] $RNA: template search produced no templates.tsv"
  exit 3
fi
echo -e "$RNA\ttemplate_search\tok\ttsv created" >> "$STATUS"

# ===== Step 3: trim + OF3 prep =====
if [[ ! -d "$TMP_BASE/templates" ]]; then
  echo -e "$RNA\tof3_template_prep\tfail\tno templates dir" >> "$STATUS"
  echo "[FAIL] $RNA: no templates/ dir from template search"
  exit 3
fi

python "$TRIM_PY" \
  --templates_tsv "$TMP_BASE/${RNA}_vfold3D_clean.mtf.templates.tsv" \
  --templates_dir "$TMP_BASE/templates" \
  --mtf "$TMP_BASE/${RNA}_vfold3D_clean.mtf" \
  --helix_pdb "$HELIX_PDB" >/dev/null 2>&1

if ! python "$PREP_PY" \
    --rna-name "$RNA" \
    --seq-file "$SEQFILE" \
    --templates-tsv "$TMP_BASE/${RNA}_vfold3D_clean.mtf.templates.tsv" \
    --templates-dir "$TMP_BASE/templates" \
    --out-dir "$TMP_BASE" ; then
  echo -e "$RNA\tof3_template_prep\tfail\tno usable motif templates" >> "$STATUS"
  echo "[FAIL] $RNA: OF3 template prep failed (no usable motif templates)"
  exit 3
fi
echo -e "$RNA\tof3_template_prep\tok\tprepared" >> "$STATUS"

# ===== Step 3.5: post-prep fixes for templates to actually work =====
# (without these, templates are silently dropped at inference time)

# 3.5a. Flatten nested CIFs: of3_structures/<m>/<m>_<chain>.cif -> of3_structures/<m>.cif
echo "[FIX] Flattening CIFs"
( cd "$TMP_BASE/of3_structures" 2>/dev/null && \
  for d in m*/; do
    [[ -d "$d" ]] || continue
    d=${d%/}
    cif=$(ls "$d"/*.cif 2>/dev/null | head -1)
    [[ -n "$cif" ]] && cp "$cif" "${d}.cif"
  done
) || { echo "[FAIL] $RNA flatten"; exit 3; }

# 3.5b. Fix .sto: ensure query has 'query_A/1-N' header (StoParser requirement)
STO="$TMP_BASE/${RNA}_templates.sto"
if [[ -f "$STO" ]] && ! grep -q "#=GS query_A/" "$STO"; then
  echo "[FIX] Adding query_A header to .sto"
  QSEQ=$(awk '/^query / {print $2; exit} /^query_A/ {print $2; exit}' "$STO")
  if [[ -n "$QSEQ" ]]; then
    QLEN=${#QSEQ}
    cp "$STO" "${STO}.before_fix"
    {
      echo "# STOCKHOLM 1.0"
      echo "#=GS query_A/1-${QLEN} DE mol:rna length:${QLEN}"
      grep "^#=GS " "${STO}.before_fix" | grep -v "query"
      awk -v qlen="$QLEN" '
        /^# STOCKHOLM/ { next }
        /^#=GS/ { next }
        /^\/\// { next }
        /^[[:space:]]*$/ { next }
        /^query / { printf "query_A/1-%d          %s\n", qlen, $2; next }
        /^query_A/ { print; next }
        { print }
      ' "${STO}.before_fix"
      echo "//"
    } > "$STO"
  fi
fi

# 3.5c. Apply v7 CIF rewriter (adds entity_poly, chem_comp, pdbx_poly_seq_scheme,
#       pdbx_audit_revision_history; renumbers residues to match cache idx_map)
echo "[FIX] Applying fix_motif_cifs_v7"
if ! python /home/dpancholi1/Projects/Openfold_template/fix_motif_cifs_v7.py \
    "$TMP_BASE/of3_structures" "$STO"; then
  echo -e "$RNA\tcif_fix\tfail\tv7 rewriter failed" >> "$STATUS"
  echo "[FAIL] $RNA: v7 CIF fix failed"
  exit 3
fi

# 3.5d. Adjust runner.yml flags (cache_directory propagation + avoid deep parse_mmcif)
RUNNER="$TMP_BASE/runner_templates_${RNA}.yml"
if [[ -f "$RUNNER" ]]; then
  sed -i 's/create_precache: false/create_precache: true/' "$RUNNER"
  sed -i 's/preparse_structures: true/preparse_structures: false/' "$RUNNER"
fi

# 3.5e. Wipe stale cache (in case this RNA was run before with broken templates)
rm -rf "$TMP_BASE/template_data/template_cache"/* 2>/dev/null
rm -rf "$TMP_BASE/template_data/template_precache"/* 2>/dev/null

echo -e "$RNA\ttemplate_fixes\tok\tCIFs/sto/yml fixed" >> "$STATUS"

# ===== Step 4: OF3 prediction with templates =====
TMP_OUT="$TMP_BASE/output_${CKBASE}"
echo "[RUN] $RNA OF3 with-template ($CKBASE) — ${NUM_DIFF}×${NUM_SEEDS}=$((NUM_DIFF*NUM_SEEDS)) preds"

if run_openfold predict \
    --query-json "$TMP_BASE/${RNA}_query_with_templates.json" \
    --inference-ckpt-path "$CKPT" \
    --runner-yaml "$TMP_BASE/runner_templates_${RNA}.yml" \
    --use-msa-server False \
    --use-templates True \
    --num-diffusion-samples "$NUM_DIFF" \
    --num-model-seeds "$NUM_SEEDS" \
    --output-dir "$TMP_OUT" ; then
  echo -e "$RNA\twithtemplate_${CKBASE}\tok\tdone" >> "$STATUS"
  echo "[DONE] $RNA OF3 prediction"
  exit 0
else
  echo -e "$RNA\twithtemplate_${CKBASE}\tfail\trun_openfold failed" >> "$STATUS"
  echo "[FAIL] $RNA: run_openfold returned non-zero"
  exit 4
fi

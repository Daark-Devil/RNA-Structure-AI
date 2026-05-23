#!/bin/bash
# Usage:
#   run_ss_msa_pipeline.sh <RNA_NAME> <SEQ_FILE> <SS_FILE> <PIPELINE_A_V2_DIR> [REF_PDB]
set -e

RNA=$1
SEQ_FILE=$2
SS_FILE=$3
SRC_DIR=$4
REF_PDB=$5

if [ -z "$RNA" ] || [ -z "$SEQ_FILE" ] || [ -z "$SS_FILE" ] || [ -z "$SRC_DIR" ]; then
    echo "Usage: $0 <RNA_NAME> <SEQ_FILE> <SS_FILE> <PIPELINE_A_V2_DIR> [REF_PDB]"
    echo
    echo "Required:"
    echo "  RNA_NAME            e.g. 3G8T"
    echo "  SEQ_FILE            path to <RNA>.seq"
    echo "  SS_FILE             path to <RNA>.2d (dot-bracket)"
    echo "  PIPELINE_A_V2_DIR   dir containing template files (e.g. .../withtemplate/3G8T)"
    echo
    echo "Optional:"
    echo "  REF_PDB             native structure for RMSD comparison"
    exit 1
fi

PIPELINE=~/ss_msa_pipeline
WORKDIR=~/ss_msa_test/$RNA
GEN=~/ss_msa_generator/ss_to_synthetic_msa.py
CKPT=/home/dpancholi1/Projects/Openfold_template/params/of3-p2-145k.pt

source ~/miniconda3/etc/profile.d/conda.sh
conda activate openfold3_rna_templ
module load cuda12.8/toolkit/12.8.1 2>/dev/null
unset CUDA_HOME
[ -x "$(command -v nvcc)" ] && export CUDA_HOME="$(dirname "$(dirname "$(command -v nvcc)")")"

echo "===================================="
echo "  SS-MSA Pipeline: $RNA"
echo "===================================="

echo
echo "[1/5] Setting up workdir at $WORKDIR"
mkdir -p $WORKDIR/synthetic_msa/chain_A
mkdir -p $WORKDIR/template_data/template_cache
mkdir -p $WORKDIR/template_data/logs
mkdir -p $WORKDIR/of3_structures

echo "[2/5] Copying template files from $SRC_DIR"
cp $SRC_DIR/${RNA}_templates.sto $WORKDIR/${RNA}_templates.sto
cp $SRC_DIR/of3_structures/*.cif $WORKDIR/of3_structures/ 2>/dev/null || true
echo "  Templates: $(ls $WORKDIR/of3_structures/*.cif 2>/dev/null | wc -l) cif files"

echo "[3/5] Generating synthetic MSA from $SS_FILE"
python $GEN \
    --seq $SEQ_FILE \
    --ss  $SS_FILE \
    --out $WORKDIR/synthetic_msa/chain_A/rfam_hits.sto \
    --query_name ${RNA}_A 2>&1 | tail -5

echo "[4/5] Building JSON and YAML configs"
SEQ=$(cat $SEQ_FILE | head -1 | tr -d '\n\r ')
TPL_IDS=$(python3 -c "
import json
d = json.load(open('$SRC_DIR/${RNA}_query_with_templates.json'))
ids = list(d['queries'].values())[0]['chains'][0]['template_entry_chain_ids']
print(json.dumps(ids))
")

sed -e "s|{{RNA}}|$RNA|g" \
    -e "s|{{SEQ}}|$SEQ|g" \
    -e "s|{{WORKDIR}}|$WORKDIR|g" \
    -e "s|{{TPL_IDS}}|$TPL_IDS|g" \
    $PIPELINE/templates/with_msa.json.tpl > $WORKDIR/${RNA}_with_synthetic_msa.json

sed -e "s|{{RNA}}|$RNA|g" \
    -e "s|{{SEQ}}|$SEQ|g" \
    -e "s|{{WORKDIR}}|$WORKDIR|g" \
    -e "s|{{TPL_IDS}}|$TPL_IDS|g" \
    $PIPELINE/templates/template_only.json.tpl > $WORKDIR/${RNA}_baseline.json

sed -e "s|{{WORKDIR}}|$WORKDIR|g" \
    $PIPELINE/templates/runner.yml.tpl > $WORKDIR/runner_${RNA}.yml

echo "  Configs created"

echo "[5/5] Running OF3 prediction (templates + synthetic MSA)"
cd $WORKDIR

run_openfold predict \
    --query-json ${RNA}_with_synthetic_msa.json \
    --inference-ckpt-path $CKPT \
    --runner-yaml runner_${RNA}.yml \
    --use-msa-server False \
    --use-templates True \
    --num-diffusion-samples 10 \
    --num-model-seeds 3 \
    --output-dir ./of3_output_msa 2>&1 | tee run_msa.log | tail -10

if [ -n "$REF_PDB" ] && [ -f "$REF_PDB" ]; then
    echo
    echo "===================================="
    echo "  RMSD vs reference: $REF_PDB"
    echo "===================================="
    KABSCH=~/of3_test_install/of3_template_deliverable/tools/compare_rna_kabsch.py
    
    RMSDS=$(for f in $WORKDIR/of3_output_msa/${RNA}_msa_test/seed_*/*sample*_model.cif; do
        python $KABSCH --ref $REF_PDB --pred $f 2>/dev/null | grep -oE '[0-9]+\.[0-9]+' | head -1
    done)
    
    BEST=$(echo "$RMSDS" | sort -g | head -1)
    MEAN=$(echo "$RMSDS" | python3 -c "import sys; v=[float(x) for x in sys.stdin if x.strip()]; print(f'{sum(v)/len(v):.3f}')" 2>/dev/null)
    
    echo "Best RMSD:  $BEST Å"
    echo "Mean RMSD:  $MEAN Å"
    echo "All values:"
    echo "$RMSDS" | sort -g
fi

echo
echo "===================================="
echo "  Pipeline complete: $RNA"
echo "===================================="
echo "Workdir: $WORKDIR"
echo "Output:  $WORKDIR/of3_output_msa/"

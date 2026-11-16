#!/usr/bin/env bash
#
# Set up and submit FragPipe MS searches for all TMT plexes.
#
# Steps:
#   1. Copy per-plex FASTAs to per_plex_fragpipe/ with rev_ decoys appended
#      (leaves per_plex/ untouched for MaxQuant, which adds decoys internally)
#   2. Generate per-plex manifests from TMT mapping
#   3. Submit SLURM array job to run FragPipe per plex
#
# Prerequisites:
#   - RAW files downloaded and organized into per-plex subdirectories:
#       RAW/<plex>/*.raw   (generate_manifests.py searches one level down)
#   - Per-plex FASTAs generated (sbatch scripts/fasta_gen/submit_proteogenomics.sh)
#   - FragPipe installed and FRAGPIPE_BIN set (or edit submit_fragpipe.sh)
#   - Workflow file exported from FragPipe GUI (TMT-10 closed search recommended)
#
# Usage:
#   bash scripts/ms_search/run_ms_search.sh [--single-fasta FASTA] [workflow_template]
#
# Arguments:
#   --single-fasta FASTA  Optional. Search every plex against this one
#                         FASTA instead of per-plex custom FASTAs. Decoys
#                         are added to a copy of this file (named
#                         {name}_decoy{ext}) and that file is used for
#                         every plex. If omitted, the default per-plex
#                         mode is used (FASTA/per_plex/{plex_id}.fasta).
#   workflow_template     Optional path to a FragPipe .workflow file.
#                      If provided, per-plex workflows are generated with
#                      the correct FASTA path patched in.
#                      If omitted, place a default workflow at:
#                        MS_SEARCH/fragpipe.workflow
#

set -euo pipefail

# --------- PATHS ----------
SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"
META_DIR="/home/lam.wayn/aas_stabilization/reprocessing/PDC_meta"
DATA_DIR="/scratch/lam.wayn"

RAW_DIR="${DATA_DIR}/RAW"                            # contains <plex>/*.raw subdirectories
FASTA_SRC_DIR="${DATA_DIR}/FASTA/per_plex"          # source (no decoys, used by MaxQuant)
FASTA_DIR="${DATA_DIR}/FASTA/per_plex_fragpipe"     # FragPipe copy (with rev_ decoys)
TMT_MAP="${META_DIR}/pdc_file_tmt_map.tsv"
SEARCH_DIR="${DATA_DIR}/MS_SEARCH"
# ---------------------------

SINGLE_FASTA_INPUT=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --single-fasta)
            SINGLE_FASTA_INPUT="$2"
            shift 2
            ;;
        *)
            break
            ;;
    esac
done
WORKFLOW_TEMPLATE="${1:-}"

module load python/3.13.5 2>/dev/null || true

echo "========================================="
echo "MS Search Pipeline Setup"
echo "========================================="
echo ""

# Validate inputs
for f in "$TMT_MAP"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required file not found: $f"
        exit 1
    fi
done

if [[ -n "$SINGLE_FASTA_INPUT" ]]; then
    if [[ ! -f "$SINGLE_FASTA_INPUT" ]]; then
        echo "ERROR: --single-fasta file not found: $SINGLE_FASTA_INPUT"
        exit 1
    fi
    if [[ ! -d "$RAW_DIR" ]]; then
        echo "ERROR: Required directory not found: $RAW_DIR"
        exit 1
    fi
else
    for d in "$RAW_DIR" "$FASTA_SRC_DIR"; do
        if [[ ! -d "$d" ]]; then
            echo "ERROR: Required directory not found: $d"
            exit 1
        fi
    done
fi

# RAW files are expected one level down: RAW_DIR/<plex>/*.raw
NUM_RAW_SUBDIR=$(find "$RAW_DIR" -mindepth 2 -maxdepth 2 -type f -iname '*.raw' | wc -l | tr -d ' ')
NUM_RAW_TOP=$(find "$RAW_DIR" -mindepth 1 -maxdepth 1 -type f -iname '*.raw' | wc -l | tr -d ' ')
if [[ "$NUM_RAW_SUBDIR" -eq 0 ]]; then
    echo "ERROR: No RAW files found in plex subdirectories of $RAW_DIR"
    echo "  Expected layout: $RAW_DIR/<plex>/*.raw"
    if [[ "$NUM_RAW_TOP" -gt 0 ]]; then
        echo "  Found $NUM_RAW_TOP RAW file(s) directly in $RAW_DIR;"
        echo "  move them into plex subdirectories first."
    fi
    exit 1
fi
echo "Found $NUM_RAW_SUBDIR RAW files in plex subdirectories of $RAW_DIR"
if [[ "$NUM_RAW_TOP" -gt 0 ]]; then
    echo "WARNING: $NUM_RAW_TOP RAW file(s) directly in $RAW_DIR will be ignored"
fi
echo ""

# Step 1: Add rev_ decoys
SINGLE_FASTA=""
if [[ -n "$SINGLE_FASTA_INPUT" ]]; then
    echo "[$(date)] Step 1: Adding decoys to single FASTA..."
    DECOY_OUT=$(python3 "${SCRIPTS_DIR}/add_decoys.py" \
        --single-fasta "$SINGLE_FASTA_INPUT" \
        --output-dir "$FASTA_DIR" | tee /dev/stderr)
    SINGLE_FASTA=$(echo "$DECOY_OUT" | sed -n 's/^OUTPUT_FASTA=//p')
    if [[ -z "$SINGLE_FASTA" ]]; then
        echo "ERROR: add_decoys.py did not report an output FASTA path"
        exit 1
    fi
else
    echo "[$(date)] Step 1: Building per_plex_fragpipe/ with decoys..."
    python3 "${SCRIPTS_DIR}/add_decoys.py" \
        --fasta-dir "$FASTA_SRC_DIR" \
        --output-dir "$FASTA_DIR"
fi
echo ""

# Step 2: Generate manifests
echo "[$(date)] Step 2: Generating FragPipe manifests..."

WORKFLOW_FLAG=""
if [[ -n "$WORKFLOW_TEMPLATE" ]]; then
    if [[ ! -f "$WORKFLOW_TEMPLATE" ]]; then
        echo "ERROR: Workflow template not found: $WORKFLOW_TEMPLATE"
        exit 1
    fi
    WORKFLOW_FLAG="--workflow-template $WORKFLOW_TEMPLATE"
fi

if [[ -n "$SINGLE_FASTA_INPUT" ]]; then
    FASTA_ARG="--single-fasta $SINGLE_FASTA"
else
    FASTA_ARG="--fasta-dir $FASTA_DIR"
fi

# python3 "${SCRIPTS_DIR}/generate_manifests.py" \
#     --tmt-map "$TMT_MAP" \
#     --raw-dir "$RAW_DIR" \
#     $FASTA_ARG \
#     $WORKFLOW_FLAG \
#     -o "$SEARCH_DIR"
# 
# echo ""

# Check plex list
PLEX_LIST="${SEARCH_DIR}/plex_list.txt"
if [[ ! -f "$PLEX_LIST" ]]; then
    echo "ERROR: Plex list not generated"
    exit 1
fi

NUM_PLEXES=$(wc -l < "$PLEX_LIST" | tr -d ' ')
if [[ "$NUM_PLEXES" -eq 0 ]]; then
    echo "No plexes ready for search. Check that:"
    echo "  1. RAW files exist in plex subdirectories: $RAW_DIR/<plex>/*.raw"
    echo "     and their names match file_name in $TMT_MAP"
    if [[ -n "$SINGLE_FASTA_INPUT" ]]; then
        echo "  2. Single FASTA is valid: $SINGLE_FASTA"
    else
        echo "  2. Per-plex FASTAs exist in $FASTA_DIR"
    fi
    exit 1
fi

# Check for workflow file
if [[ -z "$WORKFLOW_TEMPLATE" && ! -f "${SEARCH_DIR}/fragpipe.workflow" ]]; then
    echo "WARNING: No workflow file found."
    echo ""
    echo "Before submitting the search, provide a FragPipe workflow:"
    echo "  Option A: Export from FragPipe GUI -> ${SEARCH_DIR}/fragpipe.workflow"
    echo "  Option B: Re-run with template:"
    echo "    bash $0 /path/to/template.workflow"
    echo ""
    echo "Manifests have been generated. Plex list: $PLEX_LIST"
    echo "You can submit manually after providing the workflow:"
    echo "  sbatch --array=1-${NUM_PLEXES}%5 ${SCRIPTS_DIR}/submit_fragpipe.sh"
    exit 0
fi

# Step 3: Submit SLURM array job
echo "[$(date)] Submitting FragPipe jobs..."
echo "  Plexes to search: $NUM_PLEXES"
echo "  Concurrent jobs:  5"
echo ""

mkdir -p "${DATA_DIR}/logs"

if [[ -n "$SINGLE_FASTA_INPUT" ]]; then
    export SINGLE_FASTA
fi

sbatch --array=1 "${SCRIPTS_DIR}/submit_fragpipe.sh"

echo ""
echo "Monitor with: squeue -u \$USER -n fragpipe"
echo ""
echo "When complete, results will be in:"
echo "  ${SEARCH_DIR}/beneficial_mut_results/{plex_id}/"

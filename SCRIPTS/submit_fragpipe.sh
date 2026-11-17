#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --job-name=fragpipe
#SBATCH --output=/scratch/lam.wayn/logs/fragpipe_%A_%a.out
#SBATCH --error=/scratch/lam.wayn/logs/fragpipe_%A_%a.err

#
# Run FragPipe on one TMT plex per SLURM array task.
#
# Each plex is searched against its own custom FASTA database containing
# the reference proteome + plex-specific mutations + compensatory entries.
#
# Prerequisites:
#   - generate_manifests.py has been run (creates manifests/, workflows/, etc.)
#   - FragPipe, MSFragger, Philosopher, and IonQuant are installed
#

set -euo pipefail

# --------- PATHS ----------
DATA_DIR="/scratch/lam.wayn"
SEARCH_DIR="${DATA_DIR}/MS_SEARCH"
PLEX_LIST="${SEARCH_DIR}/plex_list.txt"
FASTA_DIR="${DATA_DIR}/FASTA/per_plex"
# If set (by run_ms_search.sh in single mode), use this FASTA for every
# plex instead of the per-plex FASTA in FASTA_DIR.
SINGLE_FASTA="${SINGLE_FASTA:-/scratch/lam.wayn/FASTA/per_plex_fragpipe/multi_singleton_ref_supp7_decoy.fasta}"

# FragPipe installation paths — adjust these to your installation
FRAGPIPE_BIN="${FRAGPIPE_BIN:-/home/lam.wayn/aas_stabilization/reprocessing/bin/fragpipe-24.0/bin/fragpipe}"
FRAGPIPE_TOOLS="${FRAGPIPE_TOOLS:-/home/lam.wayn/aas_stabilization/reprocessing/bin/fragpipe-24.0/tools}"
export JAVA_HOME="${JAVA_HOME:-/home/lam.wayn/aas_stabilization/reprocessing/bin/jdk-17.0.18+8}"

# Add mono to PATH for reading Thermo .raw
MONO_ENV_DIR="${MONO_ENV_DIR:-/home/lam.wayn/.conda/envs/mono_env}"
if [[ -x "${MONO_ENV_DIR}/bin/mono" ]]; then
    export PATH="${MONO_ENV_DIR}/bin:${PATH}"
fi
if ! command -v mono &>/dev/null; then
    echo "ERROR: mono not found on PATH."
    echo "  Expected it at: ${MONO_ENV_DIR}/bin/mono"
    exit 1
fi
# --------------------------

mkdir -p "${DATA_DIR}/logs"

# Get plex ID for this array task
PLEX_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PLEX_LIST")

if [[ -z "$PLEX_ID" ]]; then
    echo "ERROR: No plex at line ${SLURM_ARRAY_TASK_ID} of $PLEX_LIST"
    exit 1
fi

MANIFEST="${SEARCH_DIR}/manifests/${PLEX_ID}.fp-manifest"
if [[ -n "$SINGLE_FASTA" ]]; then
    FASTA="$SINGLE_FASTA"
else
    FASTA="${FASTA_DIR}/${PLEX_ID}.fasta"
fi
OUTDIR="${SEARCH_DIR}/beneficial_mut_results/${PLEX_ID}"
ANNOTATION="${SEARCH_DIR}/annotations/${PLEX_ID}_annotation.tsv"

# Use per-plex workflow or use default
if [[ -f "${SEARCH_DIR}/workflows/${PLEX_ID}.workflow" ]]; then
    WORKFLOW="${SEARCH_DIR}/workflows/${PLEX_ID}.workflow"
elif [[ -f "${SEARCH_DIR}/fragpipe.workflow" ]]; then
    WORKFLOW="${SEARCH_DIR}/fragpipe.workflow"
else
    echo "ERROR: No workflow file found."
    echo "  Expected: ${SEARCH_DIR}/workflows/${PLEX_ID}.workflow"
    echo "  Or:       ${SEARCH_DIR}/fragpipe.workflow"
    echo ""
    echo "  To create one:"
    echo "    1. Open FragPipe GUI, configure TMT-10 closed search"
    echo "    2. Export workflow to ${SEARCH_DIR}/fragpipe.workflow"
    echo "    3. Or re-run generate_manifests.py with --workflow-template"
    exit 1
fi

# Validate inputs
for f in "$MANIFEST" "$FASTA" "$WORKFLOW"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required file not found: $f"
        exit 1
    fi
done

# Skip if already completed (check for ratio_protein from TMT-Integrator)
if [[ -f "${OUTDIR}/tmt-report/ratio_protein_MD.tsv" ]]; then
    echo "[$(date)] Already completed: $PLEX_ID"
    exit 0
fi

mkdir -p "$OUTDIR"

echo "[$(date)] Starting FragPipe search for plex: $PLEX_ID"
echo "  Manifest:  $MANIFEST"
echo "  FASTA:     $FASTA"
echo "  Workflow:  $WORKFLOW"
echo "  Output:    $OUTDIR"
echo "  Threads:   ${SLURM_CPUS_PER_TASK}"
echo ""

# Run FragPipe in headless mode
"$FRAGPIPE_BIN" \
    --headless \
    --workflow "$WORKFLOW" \
    --manifest "$MANIFEST" \
    --workdir "$OUTDIR" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --ram "$((SLURM_MEM_PER_NODE / 1024))" \
    --config-tools-folder "$FRAGPIPE_TOOLS" \
    2>&1

echo ""
echo "[$(date)] Done: $PLEX_ID"
echo "Results: $OUTDIR"

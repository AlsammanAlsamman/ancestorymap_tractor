#!/bin/bash
#SBATCH --job-name=gwas_finemap
#SBATCH --output=logs/gwas_finemap_%j.out
#SBATCH --error=logs/gwas_finemap_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

# Create logs directory
mkdir -p logs

# Parse command line arguments
DRY_RUN=""
TARGETS=()
EXTRA_FLAGS=""
SNAKEFILE=""
CORES=""
JOBS=""
PREFER_HIGHMEM=false
HIGHMEM_MIN_MEM_MB=256000
SUBJOB_TIME="72:00:00"
SUBJOB_MIN_MEM_MB=128000

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run|-n)
            DRY_RUN="--dry-run"
            shift
            ;;
        --snakefile)
            SNAKEFILE="$2"
            shift 2
            ;;
        --cores)
            CORES="$2"
            shift 2
            ;;
        --jobs|-j)
            JOBS="$2"
            shift 2
            ;;
        --prefer-highmem)
            PREFER_HIGHMEM=true
            shift
            ;;
        --highmem-min-mem-mb)
            HIGHMEM_MIN_MEM_MB="$2"
            shift 2
            ;;
        --subjob-time)
            SUBJOB_TIME="$2"
            shift 2
            ;;
        --*)
            # Pass through other flags
            EXTRA_FLAGS="$EXTRA_FLAGS $1"
            if [[ $2 && ! $2 =~ ^-- ]]; then
                EXTRA_FLAGS="$EXTRA_FLAGS $2"
                shift 2
            else
                shift
            fi
            ;;
        *)
            # This is a target file
            TARGETS+=("$1")
            shift
            ;;
    esac
done

# Detect if this is an individual rule submission
INDIVIDUAL_RULE=false
if [[ -n "$SNAKEFILE" && "$SNAKEFILE" =~ ^rules/ ]]; then
    INDIVIDUAL_RULE=true
fi

# Set defaults based on submission type
if [[ "$INDIVIDUAL_RULE" == true ]]; then
    # Individual rule defaults
    DEFAULT_SNAKEFILE="$SNAKEFILE"
    DEFAULT_CORES="${CORES:-2}"
    DEFAULT_JOBS=1
    DEFAULT_MEM_MB=128000
    DEFAULT_MEM="128G"
    DEFAULT_TIME="72:00:00"
    DEFAULT_CPUS=2
else
    # Full pipeline defaults
    DEFAULT_SNAKEFILE="${SNAKEFILE:-Snakefile}"
    DEFAULT_CORES="${CORES:-8}"
    DEFAULT_JOBS=15
    DEFAULT_MEM_MB=128000
    DEFAULT_MEM="128G"
    DEFAULT_TIME="3-00:00:00"
    DEFAULT_CPUS=8
fi

# Default to 'all' if no targets specified and not individual rule
if [ ${#TARGETS[@]} -eq 0 ] && [[ "$INDIVIDUAL_RULE" == false ]]; then
    TARGETS=("all")
fi

# Set up environment
module load slurm
module load python/3.7.0

# Export Python path if needed
export PYTHONPATH="/s/nath-lab/alsamman/____MyCodes____/FineMappingSuite:$PYTHONPATH"

# Ensure snakemake state dirs exist before execution (helps on shared/NFS filesystems)
mkdir -p .snakemake .snakemake/log .snakemake/locks .snakemake/metadata .snakemake/incomplete

# Older Snakemake versions can race on metadata bookkeeping for clustered jobs.
# Use --drop-metadata when available to avoid false failures after successful jobs.
METADATA_FLAG=""
if snakemake --help 2>/dev/null | grep -q -- "--drop-metadata"; then
    METADATA_FLAG="--drop-metadata"
fi

# Print submission info
echo "Submission type: $([ "$INDIVIDUAL_RULE" == true ] && echo "Individual Rule" || echo "Full Pipeline")"
echo "Snakefile: $DEFAULT_SNAKEFILE"
echo "Cores: $DEFAULT_CORES"
echo "Target(s): ${TARGETS[*]:-"(from snakefile)"}"
echo "Jobs: $DEFAULT_JOBS"
echo "Sub-job time: $SUBJOB_TIME"
echo ""

# Adjust SLURM header for individual rules
if [[ "$INDIVIDUAL_RULE" == true ]]; then
    echo "Adjusting SLURM parameters for individual rule execution..."
fi

# Optional partition preference for heavy jobs: prefer highmem, fallback to serial.
PARTITION_OPT=""
if [[ "$PREFER_HIGHMEM" == true ]]; then
    if ! [[ "$HIGHMEM_MIN_MEM_MB" =~ ^[0-9]+$ ]]; then
        echo "--highmem-min-mem-mb must be a positive integer (MB): $HIGHMEM_MIN_MEM_MB" >&2
        exit 2
    fi

    if command -v sinfo >/dev/null 2>&1; then
        if sinfo -h -p highmem -t idle,mix,alloc -o "%P" 2>/dev/null | grep -q .; then
            PARTITION_OPT="--partition=highmem"
            echo "Partition preference enabled: highmem."
        else
            PARTITION_OPT="--partition=serial"
            echo "Highmem not available now; using serial partition."
        fi
    else
        PARTITION_OPT="--partition=serial"
        echo "sinfo unavailable; using serial partition."
    fi
fi

# Read cluster node restriction from analysis.yml (maintenance mode).
NODE_RESTRICT_OPT=""
if python3 -c "
import yaml, sys
cfg = yaml.safe_load(open('configs/analysis.yml'))
enabled = str(cfg.get('cluster', {}).get('restrict_nodes', {}).get('enabled', False)).lower()
sys.exit(0 if enabled == 'true' else 1)
" 2>/dev/null; then
    _partition=$(python3 -c "import yaml; cfg=yaml.safe_load(open('configs/analysis.yml')); print(cfg.get('cluster',{}).get('restrict_nodes',{}).get('partition','serial'))" 2>/dev/null)
    _nodelist=$(python3 -c "import yaml; cfg=yaml.safe_load(open('configs/analysis.yml')); print(cfg.get('cluster',{}).get('restrict_nodes',{}).get('nodelist',''))" 2>/dev/null)
    NODE_RESTRICT_OPT="--partition=${_partition:-serial}"
    if [[ -n "$_nodelist" ]]; then
        NODE_RESTRICT_OPT="$NODE_RESTRICT_OPT --nodelist=${_nodelist}"
    fi
    echo "Node restriction enabled: partition=${_partition:-serial}, nodelist=${_nodelist:-any}"
fi

CLUSTER_MEM_OPT="--mem=\$(( {resources.mem_mb} > ${SUBJOB_MIN_MEM_MB} ? {resources.mem_mb} : ${SUBJOB_MIN_MEM_MB} ))"
if [[ "$PREFER_HIGHMEM" == true ]]; then
    CLUSTER_MEM_OPT="--mem=\$(( {resources.mem_mb} > ${HIGHMEM_MIN_MEM_MB} ? {resources.mem_mb} : ${HIGHMEM_MIN_MEM_MB} ))"
    echo "Highmem memory floor active: ${HIGHMEM_MIN_MEM_MB} MB per sub-job (or rule value if higher)."
else
    echo "Sub-job memory floor active: ${SUBJOB_MIN_MEM_MB} MB per sub-job (or rule value if higher)."
fi

# Submit the workflow to SLURM
if [[ "$INDIVIDUAL_RULE" == true ]]; then
    # Individual rule submission. Honor per-rule resource blocks so long-running
    # rules are not capped by this wrapper's previous fixed 2h limit.
    snakemake \
        --cluster "sbatch \
            $PARTITION_OPT \
            $NODE_RESTRICT_OPT \
            $CLUSTER_MEM_OPT \
            --time=$SUBJOB_TIME \
            --job-name=rule_{rule}_{wildcards} \
            --cpus-per-task={resources.cores} \
            --output=logs/{rule}_{wildcards}_%j.out \
            --error=logs/{rule}_{wildcards}_%j.err" \
        --jobs ${JOBS:-$DEFAULT_JOBS} \
        --latency-wait 60 \
        --keep-going \
        --rerun-incomplete \
        --keep-incomplete \
        --default-resources mem_mb=$DEFAULT_MEM_MB time=$DEFAULT_TIME cores=$DEFAULT_CPUS tmpdir=system_tmpdir \
        --configfile configs/analysis.yml \
        --snakefile "$DEFAULT_SNAKEFILE" \
        --cores $DEFAULT_CORES \
        --verbose \
        --printshellcmds \
        --stats logs/snakemake_individual_stats.json \
        $METADATA_FLAG \
        $DRY_RUN \
        $EXTRA_FLAGS \
        "${TARGETS[@]}"
    SNAKEMAKE_EXIT_CODE=$?
else
    # Full pipeline submission - original configuration
    snakemake \
        --cluster "sbatch \
            $NODE_RESTRICT_OPT \
            $CLUSTER_MEM_OPT \
            --time=$SUBJOB_TIME \
            --job-name={rule}_{wildcards} \
            --cpus-per-task={threads} \
            --output=logs/{rule}_{wildcards}_%j.out \
            --error=logs/{rule}_{wildcards}_%j.err" \
        --jobs $DEFAULT_JOBS \
        --latency-wait 120 \
        --keep-going \
        --rerun-incomplete \
        --keep-incomplete \
        --configfile config/config.yaml \
        --snakefile "$DEFAULT_SNAKEFILE" \
        --cores $DEFAULT_CORES \
        --verbose \
        --printshellcmds \
        --stats logs/snakemake_stats.json \
        $METADATA_FLAG \
        $DRY_RUN \
        $EXTRA_FLAGS \
        "${TARGETS[@]}"
    SNAKEMAKE_EXIT_CODE=$?
fi

# Print completion message
if [[ ${SNAKEMAKE_EXIT_CODE:-1} -eq 0 ]]; then
    echo "Snakemake workflow completed for target(s): ${TARGETS[*]}"
    echo "Monitor jobs with: squeue -u $USER"
else
    echo "Snakemake workflow failed for target(s): ${TARGETS[*]}" >&2
fi

exit ${SNAKEMAKE_EXIT_CODE:-1}
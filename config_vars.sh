# Fill in custom output paths to be used for all downstream analyses

# DATA_ROOT is a folder where durable data files should be stored
DATA_ROOT="/scratch/users/bparks/BPCells_final"
# TMP_ROOT a folder where temporary data files should be stored when
#    fast disk performance is required (ideally TMP_ROOT is stored on an SSD)
TMP_ROOT="/lscratch/bparks"
# RESULTS_ROOT is a folder where timing results should be saved
RESULTS_ROOT="/oak/stanford/groups/wjg/bparks/BPCells/06_final/results"

# Give a command prefix for running via container. 
# Should contain two python-style patterns of {OMP_NUM_THREADS} and {CONTAINER}
# Containers should have read-write access to DATA_ROOT, TMP_ROOT, and RESULTS_ROOT on the host filesystem
CONTAINER_TEMPLATE="singularity run --no-home --cleanenv --env OMP_NUM_THREADS={OMP_NUM_THREADS} /oak/stanford/groups/wjg/bparks/inbox/final-containers/{CONTAINER}.sif" 


# Set to non-empty string for a laptop config
IS_LAPTOP=""


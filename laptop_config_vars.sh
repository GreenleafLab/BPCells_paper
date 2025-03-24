# Fill in custom output paths to be used for all downstream analyses

# DATA_ROOT is a folder where durable data files should be stored
DATA_ROOT="/home/bparks/Downloads/BPCells_paper/datasets"
# TMP_ROOT a folder where temporary data files should be stored when
#    fast disk performance is required (ideally TMP_ROOT is stored on an SSD)
TMP_ROOT="/tmp/BPCells_paper"
# RESULTS_ROOT is a folder where timing results should be saved
RESULTS_ROOT="/home/bparks/Downloads/BPCells_paper/results"

# Give a command prefix for running via container. 
# Should contain two python-style patterns of {OMP_NUM_THREADS} and {CONTAINER}
# Containers should have read-write access to DATA_ROOT and TMP_ROOT on the host filesystem
CONTAINER_TEMPLATE="apptainer run -B /home/bparks/Downloads/BPCells_paper --no-home --cleanenv --env OMP_NUM_THREADS={OMP_NUM_THREADS} /home/bparks/Downloads/BPCells_paper/singularity_images/{CONTAINER}.sif" 


# Set to non-empty string for a laptop config
IS_LAPTOP="true"
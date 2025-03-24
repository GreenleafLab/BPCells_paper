# BPCells paper docker setup

This folder contains a Dockerfile and software version information for creating the container image used for benchmarking.

## Downloading and running the container

For actually running the container, you'll want to download the container from Zendo where it is available in .sif and OCI tar formats: [https://doi.org/10.5281/zenodo.15066175](https://doi.org/10.5281/zenodo.15066175). The .sif file can be used directly with singularity/apptainer (see [`config_vars.sh`](../config_vars.sh) or [`laptop_config_vars.sh`](../laptop_config_vars.sh) for example usage). The tar file should be loadable via [`docker image load`](https://docs.docker.com/reference/cli/docker/image/load/). Feel free to file an issue if you have trouble with the docker setup -- the paper was run fully with singularity, so running in docker is comparatively less tested.

## Container contents

The main build script is in [`Dockerfile`](./Dockerfile). It uses [`requirements.txt`](./requirements.txt) and [`conda-env.yml`](./conda-env.yml) to install locked versions of Python dependencies, and uses a date-locked version of the Posit public package manager to try installing locked versions of R packages (though this won't lock R packages installed outside of the Posit package manager). See [`sessionInfo.txt`](./sessionInfo.txt) for specific versions installed in the container. These software versions should be current as of December 2024, though in some cases version conflicts prevent the absolute latest available version from being installed.

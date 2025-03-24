#!/bin/bash --norc

. "${MAMBA_ROOT_PREFIX}/etc/profile.d/mamba.sh"
eval "$($MAMBA_EXE shell hook --shell posix)"
micromamba activate "${ENV_NAME}"
exec "$@"
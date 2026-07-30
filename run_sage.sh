#!/usr/bin/env bash
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

if [[ -n "${DRINFELD_SAGE_PYTHON:-}" ]]; then
    exec "${DRINFELD_SAGE_PYTHON}" "$@"
fi

if command -v sage >/dev/null 2>&1; then
    exec sage -python "$@"
fi

if command -v python3 >/dev/null 2>&1 \
    && python3 -c 'import sage.all' >/dev/null 2>&1; then
    exec python3 "$@"
fi

printf '%s\n' \
    "SageMath was not found." \
    "Activate the environment from environment.yml or set DRINFELD_SAGE_PYTHON." \
    >&2
exit 127

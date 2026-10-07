#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${WM_PROJECT_DIR:-}" ]]; then
    if [[ -f /opt/openfoam12/etc/bashrc ]]; then
        set +u
        source /opt/openfoam12/etc/bashrc
        set -u
    else
        echo "ERROR: OpenFOAM-12 environment not found." >&2
        echo "Source OpenFOAM-12 first, then rerun: ./install_and_test.sh" >&2
        exit 2
    fi
fi

ROOT="$(cd "$(dirname "$0")" && pwd)"
"$ROOT/tests/static_selfcheck.sh"
"$ROOT/Allwmake"
"$ROOT/tests/smoke_cavity.sh"
if command -v mpirun >/dev/null 2>&1; then
    "$ROOT/tests/smoke_cavity_mpi.sh" "${SCHURIKEN_MPI_NP:-2}"
else
    echo "SCHURIKEN_MPI_SMOKE_TEST=SKIP reason=mpirun_not_found"
fi

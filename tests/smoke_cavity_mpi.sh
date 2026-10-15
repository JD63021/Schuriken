#!/usr/bin/env bash
# MPI smoke test: AMGsensei forms a common distributed basis and uses it.
set -euo pipefail

NP="${1:-${SCHURIKEN_MPI_NP:-2}}"
[[ "$NP" =~ ^[0-9]+$ ]] && (( NP >= 2 )) || {
    echo "ERROR: MPI smoke test requires NP >= 2" >&2
    exit 2
}

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "$ROOT/tests/lib/common.sh"
schuriken_require_openfoam
command -v mpirun >/dev/null 2>&1 || schuriken_fail "mpirun not found"
schuriken_require_library "$ROOT"

CASE="${SCHURIKEN_MPI_TEST_CASE:-/tmp/schuriken-cavity-mpi-${USER:-user}}"
schuriken_copy_cavity "$CASE"
python3 "$ROOT/tests/lib/setup_cavity.py" "$CASE" --lib "$(schuriken_library)" \
    --np "$NP"

cd "$CASE"
echo "SCHURIKEN_MPI_SMOKE_MESH=blockMesh"
blockMesh > log.blockMesh 2>&1
test -f constant/polyMesh/points || schuriken_fail "blockMesh failed"

echo "SCHURIKEN_MPI_SMOKE_DECOMPOSE=np=$NP"
decomposePar -force > log.decomposePar 2>&1
for ((r=0; r<NP; ++r)); do
    test -d "processor$r" || schuriken_fail "missing processor$r"
done

echo "SCHURIKEN_MPI_SMOKE_SOLVE=mpirun -np $NP foamRun -parallel"
mpirun -np "$NP" foamRun -parallel > log.schuriken-mpi-smoke 2>&1 \
    || schuriken_fail "mpirun failed, see $CASE/log.schuriken-mpi-smoke"

schuriken_check_log log.schuriken-mpi-smoke
grep -q "SCHURIKEN_AMGSENSEI_CONFIG.*parallel=1.*nProcs=$NP" log.schuriken-mpi-smoke \
    || schuriken_fail "MPI config marker missing"
grep -q 'SCHURIKEN_AMGSENSEI_AUGMENT' log.schuriken-mpi-smoke \
    || schuriken_fail "MPI AMGsensei did not augment"
grep 'SCHURIKEN_AMGSENSEI_CUMULATIVE' log.schuriken-mpi-smoke | tail -1

echo
echo "SCHURIKEN_MPI_SMOKE_TEST=PASS"
echo "np=$NP"
echo "case=$CASE"

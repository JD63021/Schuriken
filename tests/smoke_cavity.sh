#!/usr/bin/env bash
# Serial smoke test: pcgninja registration, AMGsensei setup and solve path.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "$ROOT/tests/lib/common.sh"
schuriken_require_openfoam
schuriken_require_library "$ROOT"

CASE="${SCHURIKEN_TEST_CASE:-/tmp/schuriken-cavity-${USER:-user}}"
schuriken_copy_cavity "$CASE"
# Keep the smoke test short; several pressure solves verify runtime
# registration, AMGsensei setup, and the PCG/GAMG solve path.
python3 "$ROOT/tests/lib/setup_cavity.py" "$CASE" --lib "$(schuriken_library)"

cd "$CASE"

echo "SCHURIKEN_SMOKE_MESH=blockMesh"
blockMesh > log.blockMesh 2>&1
test -f constant/polyMesh/points \
    || schuriken_fail "blockMesh did not create constant/polyMesh/points"

echo "SCHURIKEN_SMOKE_SOLVE=foamRun"
foamRun > log.schuriken-smoke 2>&1 \
    || schuriken_fail "foamRun failed, see $CASE/log.schuriken-smoke"
schuriken_check_log log.schuriken-smoke
grep 'SCHURIKEN_AMGSENSEI_CUMULATIVE' log.schuriken-smoke | tail -1

echo
printf 'SCHURIKEN_SMOKE_TEST=PASS\n'
printf 'case=%s\n' "$CASE"

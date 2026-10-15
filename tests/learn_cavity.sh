#!/usr/bin/env bash
# Explicit learning control: "learn no" in p and "learn yes" merged into
# pFinal must train on the pFinal solves only, the opposite of the legacy
# odd-solve rule (cavity PISO: p and pFinal alternate).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "$ROOT/tests/lib/common.sh"
schuriken_require_openfoam
schuriken_require_library "$ROOT"

CASE="${SCHURIKEN_TEST_CASE:-/tmp/schuriken-learn-${USER:-user}}"
schuriken_copy_cavity "$CASE"
python3 "$ROOT/tests/lib/setup_cavity.py" "$CASE" --lib "$(schuriken_library)" \
    --sensei learn=no --final-sensei learn=yes

cd "$CASE"
blockMesh > log.blockMesh 2>&1 || schuriken_fail "blockMesh failed"
foamRun > log.schuriken-learn 2>&1 \
    || schuriken_fail "foamRun failed, see $CASE/log.schuriken-learn"
schuriken_check_log log.schuriken-learn

grep -q 'SCHURIKEN_AMGSENSEI_CONFIG.* learn=no ' log.schuriken-learn \
    || schuriken_fail "CONFIG line does not report learn=no"

pattern="$(grep -o 'SCHURIKEN_AMGSENSEI_SOLVE.* learn=[01]' log.schuriken-learn \
    | sed 's/.*learn=//' | tr -d '\n')"
echo "learn pattern=$pattern"
[[ -n "$pattern" && "$pattern" =~ ^(01)+$ ]] \
    || schuriken_fail "expected learning on pFinal solves only (0101...)"

echo
echo "SCHURIKEN_LEARN_TEST=PASS"
echo "case=$CASE"

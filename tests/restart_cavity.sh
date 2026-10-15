#!/usr/bin/env bash
# Basis persistence: with "persist yes" the learned basis is written to the
# time directories and a restart continues with it from the first solve.
# Without persist nothing extra is written.
#   restart_cavity.sh [NP]     NP >= 2 runs decomposed with mpirun
set -euo pipefail

NP="${1:-1}"
[[ "$NP" =~ ^[0-9]+$ ]] || { echo "ERROR: NP must be a number" >&2; exit 2; }

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "$ROOT/tests/lib/common.sh"
schuriken_require_openfoam
schuriken_require_library "$ROOT"

run_case()
{
    local log="$1"
    if (( NP >= 2 )); then
        mpirun -np "$NP" foamRun -parallel > "$log" 2>&1
    else
        foamRun > "$log" 2>&1
    fi
}

basis_file()
{
    if (( NP >= 2 )); then
        echo "processor0/$1/AMGsensei_p"
    else
        echo "$1/AMGsensei_p"
    fi
}

setup()
{
    local case="$1"; shift
    schuriken_copy_cavity "$case"
    python3 "$ROOT/tests/lib/setup_cavity.py" "$case" --lib "$(schuriken_library)" \
        --np "$NP" --end-time 0.05 --write-interval 10 \
        --tolerance 1e-8 --rel-tol 0 \
        --sensei learnOddOnly=0 --sensei tailIter=2 "$@"
    (
        cd "$case"
        blockMesh > log.blockMesh 2>&1 || schuriken_fail "blockMesh failed"
        if (( NP >= 2 )); then
            decomposePar -force > log.decomposePar 2>&1 \
                || schuriken_fail "decomposePar failed"
        fi
    )
}

CASE="${SCHURIKEN_TEST_CASE:-/tmp/schuriken-restart-${USER:-user}}"

# 1. First run with persist: the basis is written at t = 0.05
setup "$CASE" --sensei persist=yes
cd "$CASE"
run_case log.run1 || schuriken_fail "first run failed, see $CASE/log.run1"
schuriken_check_log log.run1
test -f "$(basis_file 0.05)" || schuriken_fail "no basis written to $(basis_file 0.05)"
rank1="$(grep 'SCHURIKEN_AMGSENSEI_CUMULATIVE' log.run1 | tail -1 \
    | sed -n 's/.* rank=\([0-9]*\).*/\1/p')"
echo "first run: final rank=$rank1, basis file $(basis_file 0.05)"
(( rank1 >= 1 )) || schuriken_fail "first run learned no basis"

# 2. Restart from 0.05 to 0.1: the basis is read back and used immediately
sed -i -e 's/^\(startFrom[[:space:]]*\).*;/\1latestTime;/' \
       -e 's/^\(endTime[[:space:]]*\).*;/\1 0.1;/' system/controlDict
run_case log.run2 || schuriken_fail "restart failed, see $CASE/log.run2"
schuriken_check_log log.run2
grep 'SCHURIKEN_AMGSENSEI_RESTART' log.run2 || true
grep -q "SCHURIKEN_AMGSENSEI_RESTART .* rank=$rank1 .*status=PASS" log.run2 \
    || schuriken_fail "restart did not restore the rank-$rank1 basis"
grep -m1 'SCHURIKEN_AMGSENSEI_SOLVE' log.run2 | grep -q ' usedPOD=1 ' \
    || schuriken_fail "first solve after restart did not use the restored basis"
test -f "$(basis_file 0.1)" || schuriken_fail "no basis written to $(basis_file 0.1)"

# 3. Default (persist no): nothing extra is written
NOCASE="$CASE-off"
setup "$NOCASE"
cd "$NOCASE"
run_case log.run1 || schuriken_fail "persist-off run failed, see $NOCASE/log.run1"
schuriken_check_log log.run1
test -d "$(dirname "$(basis_file 0.05)")" || schuriken_fail "persist-off run wrote no time 0.05"
if [[ -e "$(basis_file 0.05)" ]]; then
    schuriken_fail "persist no: a basis file was written"
fi

echo
echo "SCHURIKEN_RESTART_TEST=PASS"
echo "np=$NP"
echo "case=$CASE"

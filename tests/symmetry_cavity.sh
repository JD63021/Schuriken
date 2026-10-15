#!/usr/bin/env bash
# Preconditioner symmetry: with a symmetric base preconditioner (DIC) the
# balanced recycled preconditioner P must be symmetric to round-off, which
# PCG requires.  AMGsensei checkSymmetry compares y'Px with x'Py for fixed
# pseudo-random x, y in every recycled solve.
#   symmetry_cavity.sh [NP]     NP >= 2 runs decomposed with mpirun
set -euo pipefail

NP="${1:-1}"
[[ "$NP" =~ ^[0-9]+$ ]] || { echo "ERROR: NP must be a number" >&2; exit 2; }

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "$ROOT/tests/lib/common.sh"
schuriken_require_openfoam
schuriken_require_library "$ROOT"

CASE="${SCHURIKEN_TEST_CASE:-/tmp/schuriken-symmetry-${USER:-user}}"
schuriken_copy_cavity "$CASE"
python3 "$ROOT/tests/lib/setup_cavity.py" "$CASE" --lib "$(schuriken_library)" \
    --np "$NP" --preconditioner DIC --end-time 0.05 --tolerance 1e-8 --rel-tol 0 \
    --sensei tailIter=2 --sensei learnOddOnly=0 --sensei checkSymmetry=yes

cd "$CASE"
LOG=log.schuriken-symmetry
blockMesh > log.blockMesh 2>&1 || schuriken_fail "blockMesh failed"
if (( NP >= 2 )); then
    command -v mpirun >/dev/null 2>&1 || schuriken_fail "mpirun not found"
    decomposePar -force > log.decomposePar 2>&1 || schuriken_fail "decomposePar failed"
    mpirun -np "$NP" foamRun -parallel > "$LOG" 2>&1 \
        || schuriken_fail "mpirun failed, see $CASE/$LOG"
else
    foamRun > "$LOG" 2>&1 || schuriken_fail "foamRun failed, see $CASE/$LOG"
fi
schuriken_check_log "$LOG"

python3 - "$LOG" <<'PY'
import re, sys
s = open(sys.argv[1]).read()
p = [float(x) for x in re.findall(r'SCHURIKEN_AMGSENSEI_SYMMETRY.*?relAsymP=(\S+)', s)]
m = [float(x) for x in re.findall(r'SCHURIKEN_AMGSENSEI_SYMMETRY.*?relAsymM=(\S+)', s)]
if not p:
    sys.exit('ERROR: no SCHURIKEN_AMGSENSEI_SYMMETRY lines')
print(f'checks={len(p)} max relAsymP={max(p):.3e} max relAsymM={max(m):.3e}')
if not max(p) < 1e-10:
    sys.exit('ERROR: recycled preconditioner is not symmetric')
PY

echo
echo "SCHURIKEN_SYMMETRY_TEST=PASS"
echo "np=$NP"
echo "case=$CASE"

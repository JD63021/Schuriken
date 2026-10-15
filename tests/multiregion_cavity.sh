#!/usr/bin/env bash
# Per-region state: foamMultiRun alternates two independent cavity regions of
# equal size but different lid speed.  Each region's p must get its own
# AMGsensei state (one RESET per region, never shared or reset again).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "$ROOT/tests/lib/common.sh"
schuriken_require_openfoam
schuriken_require_library "$ROOT"

CASE="${SCHURIKEN_TEST_CASE:-/tmp/schuriken-multiregion-${USER:-user}}"
schuriken_copy_cavity "$CASE"
python3 "$ROOT/tests/lib/setup_cavity.py" "$CASE" --lib "$(schuriken_library)" \
    --end-time 0.05 --sensei learnOddOnly=0 --sensei tailIter=2

cd "$CASE"

# Single-region files become the template of each region
fields=(0/*)
for r in cavityA cavityB; do
    mkdir -p "system/$r" "constant/$r" "0/$r"
    cp system/blockMeshDict system/fvSchemes system/fvSolution "system/$r/"
    cp constant/physicalProperties constant/momentumTransport "constant/$r/"
    cp "${fields[@]}" "0/$r/"
done
sed -i 's/uniform (1 0 0)/uniform (2 0 0)/' 0/cavityB/U
grep -q 'uniform (2 0 0)' 0/cavityB/U || schuriken_fail "could not change the cavityB lid speed"

mv system/fvSolution system/fvSolution.singleRegion
cat > system/fvSolution <<'EOF'
FoamFile
{
    format      ascii;
    class       dictionary;
    object      fvSolution;
}

PIMPLE
{
    nOuterCorrectors 1;
}
EOF

python3 - system/controlDict <<'PY'
import re, sys
from pathlib import Path
p = Path(sys.argv[1])
s = p.read_text()
s, n = re.subn(r'(?m)^solver\s+incompressibleFluid;',
               'regionSolvers\n{\n    cavityA incompressibleFluid;\n'
               '    cavityB incompressibleFluid;\n}', s)
if n != 1:
    raise SystemExit('could not find "solver incompressibleFluid;"')
s = re.sub(r'(?m)^application\s+\S+;', 'application     foamMultiRun;', s)
p.write_text(s)
PY

for r in cavityA cavityB; do
    blockMesh -region "$r" > "log.blockMesh.$r" 2>&1 \
        || schuriken_fail "blockMesh -region $r failed"
done

foamMultiRun > log.schuriken-multiregion 2>&1 \
    || schuriken_fail "foamMultiRun failed, see $CASE/log.schuriken-multiregion"
LOG=log.schuriken-multiregion
if grep -q 'FOAM FATAL' "$LOG"; then
    schuriken_fail "OpenFOAM fatal error in $LOG"
fi

grep 'SCHURIKEN_AMGSENSEI_RESET' "$LOG" || true
nReset="$(grep -c 'SCHURIKEN_AMGSENSEI_RESET' "$LOG" || true)"
nA="$(grep -c 'SCHURIKEN_AMGSENSEI_RESET region=cavityA ' "$LOG" || true)"
nB="$(grep -c 'SCHURIKEN_AMGSENSEI_RESET region=cavityB ' "$LOG" || true)"
for r in cavityA cavityB; do
    grep "SCHURIKEN_AMGSENSEI_CUMULATIVE region=$r " "$LOG" | tail -1
done

(( nReset == 2 && nA == 1 && nB == 1 )) \
    || schuriken_fail "expected one state per region, got $nReset resets (A=$nA B=$nB)"
for r in cavityA cavityB; do
    updates="$(grep "SCHURIKEN_AMGSENSEI_CUMULATIVE region=$r " "$LOG" | tail -1 \
        | sed -n 's/.* acceptedUpdates=\([0-9]*\).*/\1/p')"
    (( ${updates:-0} >= 1 )) || schuriken_fail "region $r never learned a basis"
done

echo
echo "SCHURIKEN_MULTIREGION_TEST=PASS"
echo "case=$CASE"

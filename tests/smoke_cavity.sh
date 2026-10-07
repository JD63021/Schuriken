#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${WM_PROJECT_DIR:-}" ]]; then
    echo "ERROR: OpenFOAM environment is not loaded." >&2
    echo "source /opt/openfoam12/etc/bashrc" >&2
    exit 2
fi

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

test -f "$FOAM_USER_LIBBIN/libSchuriken.so" || {
    echo "ERROR: libSchuriken.so was not installed" >&2
    echo "Run $ROOT/Allwmake first." >&2
    exit 3
}

SRC="$FOAM_TUTORIALS/incompressibleFluid/cavity"
CASE="${SCHURIKEN_TEST_CASE:-/tmp/schuriken-cavity-${USER:-user}}"
rm -rf "$CASE"
cp -a "$SRC" "$CASE"

python3 - "$CASE/system/controlDict" "$CASE/system/fvSolution" <<'PY'
from pathlib import Path
import re,sys
control=Path(sys.argv[1]); fv=Path(sys.argv[2])
c=control.read_text()
if re.search(r'(?m)^\s*libs\s*\(', c):
    raise SystemExit('smoke test expects stock cavity controlDict without libs')
insert='''\nlibs\n(\n    "libSchuriken.so"\n);\n'''
pos=c.find('// * * * *')
if pos < 0: pos=len(c)
c=c[:pos]+insert+'\n'+c[pos:]
# Keep the smoke test short; we only need several pressure solves to verify
# runtime registration, AMGsensei setup, and the PCG/GAMG solve path.
c=re.sub(r'(?m)^\s*endTime\s+[^;]+;', 'endTime         0.025;', c)
c=re.sub(r'(?m)^\s*writeInterval\s+[^;]+;', 'writeInterval   1000000;', c)
control.write_text(c)

s=fv.read_text()
pat=re.compile(r'(?ms)^\s*p\s*\{.*?^\s*\}\n\n\s*pFinal\s*\{.*?^\s*\}')
m=pat.search(s)
if not m:
    raise SystemExit('could not locate stock cavity p/pFinal blocks')
new='''    p
    {
        solver          pcgninja;
        preconditioner
        {
            preconditioner  GAMG;
            nVcycles        1;
            smoother        symGaussSeidel;
            tolerance       1e-12;
            relTol          0;
        }
        tolerance       1e-06;
        relTol          0.1;
        maxIter         500;
        AMGsensei
        {
            enabled             true;
            trainSnapshots      3;
            addRank             2;
            maxRank             8;
            recompressRank      4;
            tailIter            5;
            learnOddOnly        1;
        }
    }

    pFinal
    {
        $p;
        tolerance       1e-06;
        relTol          0;
    }'''
s=s[:m.start()]+new+s[m.end():]
fv.write_text(s)
PY

cd "$CASE"

echo "SCHURIKEN_SMOKE_MESH=blockMesh"
blockMesh 2>&1 | tee log.blockMesh

test -f constant/polyMesh/points || {
    echo "ERROR: blockMesh did not create constant/polyMesh/points" >&2
    exit 4
}

echo "SCHURIKEN_SMOKE_SOLVE=foamRun"
foamRun 2>&1 | tee log.schuriken-smoke

grep -q 'pcgninja:  Solving for p' log.schuriken-smoke || {
    echo "ERROR: pcgninja did not appear in solver output" >&2
    exit 5
}
grep -q 'SCHURIKEN_AMGSENSEI_CONFIG' log.schuriken-smoke || {
    echo "ERROR: AMGsensei configuration marker missing" >&2
    exit 6
}
if grep -q 'FOAM FATAL' log.schuriken-smoke; then
    echo "ERROR: OpenFOAM fatal error in smoke run" >&2
    exit 7
fi

echo
printf 'SCHURIKEN_SMOKE_TEST=PASS\n'
printf 'case=%s\n' "$CASE"

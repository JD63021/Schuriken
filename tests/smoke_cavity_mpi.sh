#!/usr/bin/env bash
set -euo pipefail

NP="${1:-${SCHURIKEN_MPI_NP:-2}}"
[[ "$NP" =~ ^[0-9]+$ ]] && (( NP >= 2 )) || {
    echo "ERROR: MPI smoke test requires NP >= 2" >&2
    exit 2
}

if [[ -z "${WM_PROJECT_DIR:-}" ]]; then
    echo "ERROR: OpenFOAM environment is not loaded." >&2
    echo "source /opt/openfoam12/etc/bashrc" >&2
    exit 3
fi
command -v mpirun >/dev/null 2>&1 || { echo "ERROR: mpirun not found" >&2; exit 4; }

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
test -f "$FOAM_USER_LIBBIN/libSchuriken.so" || {
    echo "ERROR: libSchuriken.so was not installed" >&2
    echo "Run $ROOT/Allwmake first." >&2
    exit 5
}

SRC="$FOAM_TUTORIALS/incompressibleFluid/cavity"
CASE="${SCHURIKEN_MPI_TEST_CASE:-/tmp/schuriken-cavity-mpi-${USER:-user}}"
rm -rf "$CASE"
cp -a "$SRC" "$CASE"

python3 - "$CASE/system/controlDict" "$CASE/system/fvSolution" "$NP" <<'PY_MPI'
from pathlib import Path
import re,sys
control=Path(sys.argv[1]); fv=Path(sys.argv[2]); np=int(sys.argv[3])
c=control.read_text()
if not re.search(r'(?m)^\s*libs\s*\(', c):
    insert='\nlibs\n(\n    "libSchuriken.so"\n);\n'
    pos=c.find('// * * * *')
    if pos < 0: pos=len(c)
    c=c[:pos]+insert+'\n'+c[pos:]
c=re.sub(r'(?m)^\s*endTime\s+[^;]+;', 'endTime         0.025;', c)
c=re.sub(r'(?m)^\s*writeInterval\s+[^;]+;', 'writeInterval   1000000;', c)
control.write_text(c)

s=fv.read_text()
pat=re.compile(r'(?ms)^\s*p\s*\{.*?^\s*\}\n\n\s*pFinal\s*\{.*?^\s*\}')
m=pat.search(s)
if not m: raise SystemExit('could not locate stock cavity p/pFinal blocks')
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

Path(fv.parent/'decomposeParDict').write_text(f'''FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "system";
    object      decomposeParDict;
}}
numberOfSubdomains {np};
method          simple;
simpleCoeffs
{{
    n           ({np} 1 1);
    delta       0.001;
}}
''')
PY_MPI

cd "$CASE"
echo "SCHURIKEN_MPI_SMOKE_MESH=blockMesh"
blockMesh > log.blockMesh 2>&1
test -f constant/polyMesh/points || { echo "ERROR: blockMesh failed" >&2; exit 6; }

echo "SCHURIKEN_MPI_SMOKE_DECOMPOSE=np=$NP"
decomposePar -force > log.decomposePar 2>&1
for ((r=0; r<NP; ++r)); do test -d "processor$r" || { echo "ERROR: missing processor$r" >&2; exit 7; }; done

echo "SCHURIKEN_MPI_SMOKE_SOLVE=mpirun -np $NP foamRun -parallel"
mpirun -np "$NP" foamRun -parallel 2>&1 | tee log.schuriken-mpi-smoke

grep -q 'pcgninja:  Solving for p' log.schuriken-mpi-smoke || { echo "ERROR: pcgninja missing" >&2; exit 8; }
grep -q "SCHURIKEN_AMGSENSEI_CONFIG.*parallel=1.*nProcs=$NP" log.schuriken-mpi-smoke || { echo "ERROR: MPI config marker missing" >&2; exit 9; }
grep -q 'SCHURIKEN_AMGSENSEI_AUGMENT' log.schuriken-mpi-smoke || { echo "ERROR: MPI AMGsensei did not augment" >&2; exit 10; }
if grep -q 'FOAM FATAL' log.schuriken-mpi-smoke; then echo "ERROR: FOAM FATAL in MPI smoke" >&2; exit 11; fi

echo
echo "SCHURIKEN_MPI_SMOKE_TEST=PASS"
echo "np=$NP"
echo "case=$CASE"

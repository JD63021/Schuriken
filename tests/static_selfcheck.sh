#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

bash -n Allwmake Allclean install_and_test.sh tests/lib/common.sh tests/*.sh
python3 -c 'import ast, sys; ast.parse(open(sys.argv[1]).read())' tests/lib/setup_cavity.py

grep -q 'TypeName("pcgninja")' src/solver/PCGNinja.H
grep -q 'addsymMatrixConstructorToTable<PCGNinja>' src/solver/PCGNinja.C
grep -q 'updateSenseiAfterSolve' src/sensei/AMGsenseiTrainer.C
grep -q 'BATCHED_GLOBAL_POD' src/solver/PCGNinja.C
grep -q 'MPI_Allreduce' src/linalg/AllReduce.C
grep -q 'listCombineGather' src/linalg/AllReduce.C
grep -q 'class RecycledGAMG' src/preconditioner/RecycledGAMG.H
grep -q 'solver[[:space:]]*pcgninja' examples/fvSolution.pcgninja
grep -q 'libSchuriken.so' examples/controlDict.libs

# Guard against collapsing the implementation back into one giant solver file.
lines="$(wc -l < src/solver/PCGNinja.C)"
if (( lines > 450 )); then
    echo "ERROR: PCGNinja.C grew to $lines lines; keep trainer/linalg/preconditioner modular."
    exit 2
fi

for f in \
  src/solver/PCGNinja.C \
  src/sensei/AMGsensei.C \
  src/sensei/AMGsenseiState.C \
  src/sensei/AMGsenseiTrainer.C \
  src/preconditioner/RecycledGAMG.C \
  src/linalg/AllReduce.C \
  src/linalg/DistributedOps.C \
  src/linalg/DenseOps.C \
  src/linalg/PODBasis.C \
  src/linalg/ReducedSpace.C; do
    grep -q "${f#src/}" src/Make/files
  done

echo "SCHURIKEN_STATIC_SELFCHECK=PASS"
echo "pcgninja_lines=$lines"

#!/usr/bin/env bash
# Forced recompression: a small maxRank makes AMGsensei reach its rank limit
# early in a short cavity run.  The GAMG-informed spectral recompression must
# then run (not be skipped) and shrink the basis to recompressRank.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "$ROOT/tests/lib/common.sh"
schuriken_require_openfoam
schuriken_require_library "$ROOT"

CASE="${SCHURIKEN_TEST_CASE:-/tmp/schuriken-recompress-${USER:-user}}"
schuriken_copy_cavity "$CASE"
python3 "$ROOT/tests/lib/setup_cavity.py" "$CASE" --lib "$(schuriken_library)" \
    --end-time 0.1 --tolerance 1e-8 --rel-tol 0 \
    --sensei trainSnapshots=2 --sensei addRank=2 \
    --sensei maxRank=4 --sensei recompressRank=2 --sensei recompress=yes \
    --sensei tailIter=2 --sensei learnOddOnly=0

cd "$CASE"
blockMesh > log.blockMesh 2>&1 || schuriken_fail "blockMesh failed"
foamRun > log.schuriken-recompress 2>&1 \
    || schuriken_fail "foamRun failed, see $CASE/log.schuriken-recompress"
schuriken_check_log log.schuriken-recompress

LOG=log.schuriken-recompress
nDone="$(grep -c 'SCHURIKEN_AMGSENSEI_RECOMPRESS solve=.*method=amg_spectral' "$LOG" || true)"
nSkip="$(grep -c 'SCHURIKEN_AMGSENSEI_RECOMPRESS_SKIP\|method=keep_newest' "$LOG" || true)"
grep -m3 'SCHURIKEN_AMGSENSEI_AMG_SPECTRUM\|SCHURIKEN_AMGSENSEI_RECOMPRESS' "$LOG" || true
echo "recompressions=$nDone skipped_or_fallback=$nSkip"

(( nDone >= 1 )) || schuriken_fail "spectral recompression never ran"
(( nSkip == 0 )) || schuriken_fail "spectral recompression was skipped or fell back $nSkip times"
if grep 'SCHURIKEN_AMGSENSEI_RECOMPRESS solve=' "$LOG" | grep -vq 'rankNew=2 '; then
    schuriken_fail "recompression did not shrink the basis to recompressRank=2"
fi

# Default (recompress no): the full basis is kept and learning stops
cd ..
NOCASE="$CASE-off"
schuriken_copy_cavity "$NOCASE"
python3 "$ROOT/tests/lib/setup_cavity.py" "$NOCASE" --lib "$(schuriken_library)" \
    --end-time 0.1 --tolerance 1e-8 --rel-tol 0 \
    --sensei trainSnapshots=2 --sensei addRank=2 \
    --sensei maxRank=4 --sensei recompressRank=2 \
    --sensei tailIter=2 --sensei learnOddOnly=0
cd "$NOCASE"
blockMesh > log.blockMesh 2>&1 || schuriken_fail "blockMesh failed"
foamRun > "$LOG" 2>&1 || schuriken_fail "foamRun failed, see $NOCASE/$LOG"
schuriken_check_log "$LOG"
grep -q 'SCHURIKEN_AMGSENSEI_BASIS_FULL.*rank=4 ' "$LOG" \
    || schuriken_fail "recompress no: basis-full marker missing"
if grep -q 'SCHURIKEN_AMGSENSEI_RECOMPRESS' "$LOG"; then
    schuriken_fail "recompress no: recompression ran"
fi
grep 'SCHURIKEN_AMGSENSEI_CUMULATIVE' "$LOG" | tail -1 | grep -q ' rank=4 ' \
    || schuriken_fail "recompress no: final rank is not maxRank=4"

echo
echo "SCHURIKEN_RECOMPRESS_TEST=PASS"
echo "case=$CASE"

#!/usr/bin/env bash
# Build and run the serial dense/distributed linear-algebra unit tests.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "$ROOT/tests/lib/common.sh"
schuriken_require_openfoam
schuriken_require_library "$ROOT"

cd "$ROOT/tests/unit"
wmake > "$ROOT/tests/unit/log.wmake" 2>&1 \
    || schuriken_fail "unit test build failed, see tests/unit/log.wmake"
"$FOAM_USER_APPBIN/Test-schurikenLinalg"

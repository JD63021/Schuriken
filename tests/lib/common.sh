# Schuriken test helpers, sourced by the scripts in tests/.

schuriken_fail()
{
    echo "ERROR: $*" >&2
    exit 1
}

schuriken_require_openfoam()
{
    if [[ -z "${WM_PROJECT_DIR:-}" ]]; then
        echo "ERROR: OpenFOAM environment is not loaded." >&2
        echo "source /opt/openfoam12/etc/bashrc" >&2
        exit 2
    fi
}

# Library loaded by the test cases; SCHURIKEN_LIB selects another build.
schuriken_library()
{
    echo "${SCHURIKEN_LIB:-libSchuriken.so}"
}

schuriken_require_library()
{
    local lib
    lib="$(schuriken_library)"
    test -f "$FOAM_USER_LIBBIN/$lib" || {
        echo "ERROR: $lib was not installed in $FOAM_USER_LIBBIN" >&2
        echo "Run $1/Allwmake first." >&2
        exit 3
    }
}

# Copy the stock cavity tutorial to $1.  A previous case at that path is
# moved aside to <case>.oldN, never deleted.
schuriken_copy_cavity()
{
    local case="$1"
    if [[ -e "$case" ]]; then
        local n=1
        while [[ -e "$case.old$n" ]]; do n=$((n+1)); done
        mv "$case" "$case.old$n"
        echo "SCHURIKEN_TEST_PREVIOUS_CASE=$case.old$n"
    fi
    cp -a "$FOAM_TUTORIALS/incompressibleFluid/cavity" "$case"
}

schuriken_check_log()
{
    local log="$1"
    grep -q 'pcgninja:  Solving for p' "$log" \
        || schuriken_fail "pcgninja did not appear in $log"
    grep -q 'SCHURIKEN_AMGSENSEI_CONFIG' "$log" \
        || schuriken_fail "AMGsensei configuration marker missing in $log"
    if grep -q 'FOAM FATAL' "$log"; then
        schuriken_fail "OpenFOAM fatal error in $log"
    fi
}

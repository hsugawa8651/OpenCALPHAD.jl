#!/bin/bash
#
# Generate reference data using openCALPHAD
#
# Prerequisites:
#   - openCALPHAD (oc6) installed and in PATH
#   - Run from reftest/ directory
#
# Usage:
#   ./generate_reference.sh <test_name>
#   ./generate_reference.sh agcu_equil_1000K
#   ./generate_reference.sh all
#

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
MACRO_DIR="$SCRIPT_DIR/macros"
REF_DIR="$SCRIPT_DIR/reference"

# Check oc6 availability
check_oc6() {
    if ! command -v oc6 &> /dev/null; then
        echo "Error: oc6 (openCALPHAD) not found in PATH"
        echo ""
        echo "Installation instructions:"
        echo "  1. Clone: git clone https://github.com/sundmanbo/opencalphad"
        echo "  2. Build: cd opencalphad && make"
        echo "  3. Add to PATH: export PATH=\$PATH:/path/to/opencalphad/bin"
        exit 1
    fi

    echo "Found oc6: $(which oc6)"
    oc6 --version 2>/dev/null || echo "(version check not available)"
}

# Run a single test
run_test() {
    local name="$1"
    local macro="$MACRO_DIR/${name}.OCM"
    local log="$REF_DIR/${name}.log"

    if [ ! -f "$macro" ]; then
        echo "Error: Macro file not found: $macro"
        return 1
    fi

    echo "Running: $name"
    echo "  Macro: $macro"
    echo "  Output: $log"

    cd "$SCRIPT_DIR"
    oc6 < "$macro" > "$log" 2>&1

    echo "  Status: Done"
    echo ""
    echo "Next steps:"
    echo "  1. Review output: cat $log"
    echo "  2. Extract values and update: reference/${name}.json"
    echo "  3. Use parse_log.jl to assist: julia parse_log.jl $log"
}

# List available tests
list_tests() {
    echo "Available reference tests:"
    for macro in "$MACRO_DIR"/*.OCM; do
        if [ -f "$macro" ]; then
            name=$(basename "$macro" .OCM)
            echo "  - $name"
        fi
    done
}

# Main
case "${1:-}" in
    "")
        echo "Usage: $0 <test_name|all|list>"
        echo ""
        list_tests
        ;;
    "list")
        list_tests
        ;;
    "all")
        check_oc6
        for macro in "$MACRO_DIR"/*.OCM; do
            if [ -f "$macro" ]; then
                name=$(basename "$macro" .OCM)
                run_test "$name"
            fi
        done
        ;;
    *)
        check_oc6
        run_test "$1"
        ;;
esac

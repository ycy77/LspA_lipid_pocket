#!/usr/bin/env bash
set -euo pipefail
LSP_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$LSP_ROOT/config/paths.sh"
if [[ $# -eq 0 || "$1" == --help ]]; then
    echo 'Usage: bash run.sh <stage/script.py|stage/script.sh> [arguments]'
    echo 'Example: bash run.sh design/generate_backbones.py --help'
    find "$LSP_ROOT/scripts" -type f \( -name '*.py' -o -name '*.sh' \) | sort
    exit 0
fi
entry="$1"
shift
[[ "$entry" != /* && "$entry" != *..* ]] || { echo 'Use a path under scripts/' >&2; exit 2; }
script="$LSP_ROOT/scripts/$entry"
[[ -f "$script" ]] || { echo "Unknown script: $entry" >&2; exit 2; }
case "$script" in
    *.py) exec "${LSP_PYTHON:-python3}" "$script" "$@" ;;
    *.sh) exec bash "$script" "$@" ;;
esac

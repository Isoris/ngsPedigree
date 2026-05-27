#!/usr/bin/env bash
# DROP.sh — land the HPP / relatedness_atlas registry into atlas-core.
#
# Usage:
#   bash docs/hpp/atlas_core_registry/DROP.sh <ATLAS_CORE_ROOT>
#
# where ATLAS_CORE_ROOT is the local checkout of the atlas-core repo
# (the one that contains toolkit_registries/). The tarball's internal
# layout is toolkit_registries/relatedness/01_registry/*.jsonl, so it
# extracts directly into the expected location.
#
# Idempotent: re-runs overwrite the four JSONL files in place. No other
# atlas-core files are touched.

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $(basename "$0") <ATLAS_CORE_ROOT>" >&2
    exit 64
fi

ATLAS_CORE_ROOT="$1"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARBALL="${HERE%/atlas_core_registry}/atlas_core_registry.tar.gz"

if [[ ! -d "$ATLAS_CORE_ROOT" ]]; then
    echo "atlas-core root not found: $ATLAS_CORE_ROOT" >&2
    exit 66
fi
if [[ ! -d "$ATLAS_CORE_ROOT/toolkit_registries" ]]; then
    echo "warning: $ATLAS_CORE_ROOT does not contain toolkit_registries/" >&2
    echo "         continuing anyway; tarball will create it" >&2
fi
if [[ ! -f "$TARBALL" ]]; then
    echo "tarball not found: $TARBALL" >&2
    exit 66
fi

echo "Extracting $TARBALL → $ATLAS_CORE_ROOT/"
tar -xzvf "$TARBALL" -C "$ATLAS_CORE_ROOT"

TARGET_DIR="$ATLAS_CORE_ROOT/toolkit_registries/relatedness/01_registry"
echo
echo "Landed files:"
ls -la "$TARGET_DIR"

echo
echo "Running local smoke test ..."
python3 - "$TARGET_DIR" <<'PY'
import json, sys
from pathlib import Path

d = Path(sys.argv[1])
load = lambda n: [json.loads(l) for l in (d/n).read_text().splitlines() if l.strip()]

modules  = load("module_registry.jsonl")
analyses = load("analysis_registry.jsonl")
modes    = load("analysis_modes.jsonl")
layers   = load("layer_registry.jsonl")

module_names  = {m["module_name"] for m in modules}
analysis_ids  = {a["analysis_id"] for a in analyses}
analysis_prod = {a["analysis_id"]: set(a["produces"]) for a in analyses}
layer_ids     = {l["layer_id"] for l in layers}

ok = True
for m in modes:
    if m["analysis_type"] not in analysis_ids:
        print(f"FAIL: mode {m['analysis_type']} not in analysis_registry"); ok=False
    if m["produces"] not in analysis_prod.get(m["analysis_type"], set()):
        print(f"FAIL: mode {m['analysis_type']} produces {m['produces']} not declared"); ok=False
    if m["module_name"] not in module_names:
        print(f"FAIL: mode {m['analysis_type']} module {m['module_name']} not in module_registry"); ok=False
for lid in {p for a in analyses for p in a["produces"]}:
    if lid not in layer_ids:
        print(f"WARN: declared layer {lid} not in layer_registry"); ok=False

print(f"\nrows: {len(modules)} modules / {len(analyses)} analyses / {len(modes)} modes / {len(layers)} layers")
print("SMOKE TEST: PASS" if ok else "SMOKE TEST: FAIL")
sys.exit(0 if ok else 1)
PY

echo
echo "DROP complete. Commit the JSONL files in $ATLAS_CORE_ROOT as appropriate."

#!/bin/bash
# Adds MSExperiment.h to files that got it transitively through ConsensusMap.h.
# ConsensusMap.h previously included MSExperiment.h — now it only forward-declares it.

OPENMS_SRC="/home/blueicecream/openms-development/OpenMS/src"
OPENMS_ROOT="/home/blueicecream/openms-development/OpenMS"

add_include_after_consensusmap() {
    local file="$1"
    local new_include="$2"
    if grep -q "$new_include" "$file"; then
        return
    fi
    # Insert after the ConsensusMap.h include line
    sed -i "/#include.*ConsensusMap\.h/a #include $new_include" "$file"
    echo "  + added $new_include to $(basename $file)"
}

FILES=$(grep -rl "ConsensusMap\.h" "$OPENMS_SRC" "$OPENMS_ROOT/doc" --include="*.cpp" 2>/dev/null)

echo "Processing $(echo "$FILES" | wc -l) files..."

for f in $FILES; do
    needs_msexperiment=false
    grep -q "PeakMap\|MSExperiment" "$f" && needs_msexperiment=true

    if [ "$needs_msexperiment" = true ]; then
        add_include_after_consensusmap "$f" "<OpenMS/KERNEL/MSExperiment.h>"
    fi
done

echo "Done."

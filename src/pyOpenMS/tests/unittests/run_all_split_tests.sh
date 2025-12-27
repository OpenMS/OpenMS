#!/bin/bash
# Test runner for all split test files
# Part of Issue #8567 - test000.py splitting

cd "$(dirname "$0")"
source ../../../../venv_test/bin/activate

echo "======================================================================="
echo "Running All Split Test Files"
echo "======================================================================="
echo ""

total_files=0
passed_files=0
failed_files=0

test_files=(
    "test_chemistry.py"
    "test_core_structures.py"
    "test_spectrum_experiment.py"
    "test_features.py"
    "test_file_io.py"
    "test_identification.py"
    "test_targeted_mrm.py"
    "test_algorithms.py"
    "test_metadata_instrument.py"
    "test_quantification.py"
    "test_miscellaneous.py"
    "test_validation_signal.py"
)

for test_file in "${test_files[@]}"; do
    total_files=$((total_files + 1))
    echo "=== Testing: $test_file ==="
    
    if timeout 30 python "$test_file" > "/tmp/${test_file%.py}_output.txt" 2>&1; then
        # Count test results
        passed=$(grep -c "✅" "/tmp/${test_file%.py}_output.txt" || echo 0)
        failed=$(grep -c "❌" "/tmp/${test_file%.py}_output.txt" || echo 0)
        
        echo "  ✅ File completed: $passed passed, $failed failed"
        passed_files=$((passed_files + 1))
    else
        exit_code=$?
        if [ $exit_code -eq 124 ]; then
            echo "  ⚠️  File timed out (>30s)"
        else
            echo "  ❌ File failed with exit code: $exit_code"
        fi
        failed_files=$((failed_files + 1))
        
        # Show last 10 lines of output
        echo "  Last output:"
        tail -10 "/tmp/${test_file%.py}_output.txt" | sed 's/^/    /'
    fi
    echo ""
done

echo "======================================================================="
echo "Summary"
echo "======================================================================="
echo "Total files: $total_files"
echo "Completed: $passed_files"
echo "Failed/Timeout: $failed_files"
echo ""

if [ $failed_files -eq 0 ]; then
    echo "✅ All test files completed successfully!"
    exit 0
else
    echo "⚠️  Some test files had issues"
    exit 1
fi

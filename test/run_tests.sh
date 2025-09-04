#!/bin/bash

# Kraken KMC Integration Test Suite Runner
# This script must be run from the project root directory as: test/run_tests.sh

set -e  # Exit on any error

# Check if we're being run from the project root directory
if [ ! -f "src/kmc_to_jellyfish" ]; then
    echo "ERROR: This script must be run from the project root directory"
    echo ""
    echo "Usage: test/run_tests.sh"
    echo ""
    echo "Current directory: $(pwd)"
    echo "Expected to find: src/kmc_to_jellyfish"
    echo ""
    echo "Please run this script from the project root directory."
    exit 1
fi

echo "Running Kraken KMC Integration Tests..."
echo "======================================"

# Run the test suite
python test/run_tests.py

echo ""
echo "Tests completed!"

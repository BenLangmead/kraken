#!/usr/bin/env python3

"""
Test driver script for Kraken KMC integration tests.

This script runs all the test suites to validate the KMC to Jellyfish conversion
and integration with the Kraken database building workflow.
"""

import os
import sys
import subprocess
import time
from pathlib import Path

def run_test(test_name: str, test_script: str) -> bool:
    """Run a single test and return True if it passes."""
    print(f"\n{'='*60}")
    print(f"Running {test_name}")
    print(f"{'='*60}")
    
    start_time = time.time()
    
    try:
        # Run the test script
        result = subprocess.run([sys.executable, test_script], 
                              cwd=Path(__file__).parent,
                              capture_output=False)
        
        test_time = time.time() - start_time
        
        if result.returncode == 0:
            print(f"\n{test_name} PASSED ({test_time:.2f}s)")
            return True
        else:
            print(f"\n{test_name} FAILED ({test_time:.2f}s)")
            return False
            
    except Exception as e:
        test_time = time.time() - start_time
        print(f"\n{test_name} ERROR ({test_time:.2f}s): {e}")
        return False

def main():
    """Run all tests and report results."""
    print("Kraken KMC Integration Test Suite")
    print("=" * 60)
    
    # Get the test directory
    test_dir = Path(__file__).parent
    
    # Define all tests to run
    tests = [
        ("KMC Conversion Tool Tests", "test_kmc_conversion.py"),
        ("Integration Tests", "test_integration.py"),
    ]
    
    # Check if all test files exist
    missing_tests = []
    for test_name, test_script in tests:
        if not (test_dir / test_script).exists():
            missing_tests.append(test_script)
    
    if missing_tests:
        print(f"ERROR: Missing test files: {missing_tests}")
        sys.exit(1)
    
    # Run all tests
    passed = 0
    failed = 0
    
    for test_name, test_script in tests:
        if run_test(test_name, test_script):
            passed += 1
        else:
            failed += 1
    
    # Report final results
    print(f"\n{'='*60}")
    print("TEST SUITE RESULTS")
    print(f"{'='*60}")
    print(f"Total tests: {len(tests)}")
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")
    print(f"Success rate: {passed/(passed+failed)*100:.1f}%")
    
    if failed == 0:
        print("\nAll tests PASSED!")
        sys.exit(0)
    else:
        print(f"\n{failed} test(s) FAILED")
        sys.exit(1)

if __name__ == "__main__":
    main()

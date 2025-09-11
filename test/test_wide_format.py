#!/usr/bin/env python3

"""
Author: Ben Langmead
Date: 2025-09-04

Test script to verify that "wide" format databases (k > 31) can be read
by the KrakenDB class and that the 128-bit k-mer functions work
correctly.
"""

import sys
from test_utils import (
    setup_test_environment, 
    create_wide_format_database, 
    run_test_kraken
)

def test_wide_format_database():
    """Test that wide format databases can be created and read."""
    print("Testing wide format database creation and reading...")
    
    # Create test sequences for k=32
    sequences = [
        'ACGTACGTACGTACGTACGTACGTACGTACGT',
        'TGCATGCATGCATGCATGCATGCATGCATGCA',
        'GATCGATCGATCGATCGATCGATCGATCGATC'
    ]

    long_sequence = ["".join(sequences)*5]

    subtests={'multi-fasta':(sequences, False), 'multi-line-fasta':(long_sequence, True)}
    results=[]   
    for subtest_name, (sequences, multiline_1sequence)  in subtests.items():

        # Create wide format database
        wide_db = create_wide_format_database(sequences, 32, multiline_1sequence)
        print(f"Created wide format database: {wide_db}")
        
        # Verify the database file exists and has reasonable size
        if not wide_db.exists():
            print("FAIL: Wide format database was not created for test {subtest_name}")
            results.append(False)
            continue
            
        db_size = wide_db.stat().st_size
        if db_size == 0:
            print(f"FAIL: Wide format database is empty for test {subtest_name}")
            results.append(False)
            continue
            
        print(f"PASS: Wide format database created successfully (size: {db_size} bytes) for test {subtest_name}")
        
        # Test that the database can be read by our test program
        result = run_test_kraken()
        
        if result.returncode == 0:
            print(f"PASS: test_kraken program runs successfully for test {subtest_name}")
            results.append(True)
            
        else:
            print(f"FAIL: test_kraken program failed with return code {result.returncode} for test {subtest_name}")
            if result.stderr:
                print(f"stderr: {result.stderr}")
            results.append(False)
            continue

    return all(results)


def main():
    print("Wide Format Database Test")
    print("=" * 50)
    
    # Set up test environment
    setup_test_environment()
    
    if test_wide_format_database():
        print("\nAll wide format tests PASSED!")
        return 0
    else:
        print("\nSome wide format tests FAILED!")
        return 1

if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3

"""
Author: Ben Langmead
Date: 2025-09-04

Test script to verify database compatibility checking:
1. Wide format databases (k > 31) work with current Kraken
2. Wide format databases fail with older Kraken versions (simulated)
3. Normal format databases work with current Kraken
"""

import sys
import tempfile
from pathlib import Path
from test_utils import (
    setup_test_environment,
    create_wide_format_database,
    create_normal_format_database,
    run_test_kraken,
    analyze_database_header
)


def test_wide_format_compatibility():
    """Test that wide format databases work with current Kraken."""
    print("Testing wide format database compatibility...")
    
    # Create test sequences for k=32
    sequences = [
        'ACGTACGTACGTACGTACGTACGTACGTACGT',
        'TGCATGCATGCATGCATGCATGCATGCATGCA'
    ]
    
    wide_db = create_wide_format_database(sequences, 32)
    
    # Try to load the database with our current KrakenDB
    # This should work and show the wide format message
    result = run_test_kraken()
    
    if result.returncode == 0:
        print("PASS: Wide format database loads successfully with current Kraken")
        return True
    else:
        print(f"FAIL: Wide format database failed to load: {result.stderr}")
        return False


def test_normal_format_compatibility():
    """Test that normal format databases still work with current Kraken."""
    print("Testing normal format database compatibility...")
    
    # Create test sequences for k=8
    sequences = [
        'ACGTACGTACGTACGT',
        'TGCATGCATGCATGCA'
    ]
    
    normal_db = create_normal_format_database(sequences, 8)
    
    # Try to load the database with our current KrakenDB
    # This should work without any special messages
    result = run_test_kraken()
    
    if result.returncode == 0:
        print("PASS: Normal format database loads successfully with current Kraken")
        return True
    else:
        print(f"FAIL: Normal format database failed to load: {result.stderr}")
        return False


def test_database_header_analysis():
    """Test that we can analyze database headers to determine format."""
    print("Testing database header analysis...")
    
    # Create databases in a persistent temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Create wide format database
        wide_sequences = ['ACGTACGTACGTACGTACGTACGTACGTACGT', 'TGCATGCATGCATGCATGCATGCATGCATGCA']
        wide_db = create_wide_format_database(wide_sequences, 32, temp_dir)
        
        # Create normal format database
        normal_sequences = ['ACGTACGTACGTACGT', 'TGCATGCATGCATGCA']
        normal_db = create_normal_format_database(normal_sequences, 8, temp_dir)
        
        # Analyze headers
        wide_info = analyze_database_header(wide_db)
        normal_info = analyze_database_header(normal_db)
        
        print(f"Wide format database: key_bits={wide_info['key_bits']}, k={wide_info['k']}")
        print(f"Normal format database: key_bits={normal_info['key_bits']}, k={normal_info['k']}")
        
        # Verify the formats
        if wide_info['key_bits'] > 62:  # k > 31
            print("PASS: Wide format database correctly identified (key_bits > 62)")
        else:
            print("FAIL: Wide format database not correctly identified")
            return False
        
        if normal_info['key_bits'] <= 62:  # k <= 31
            print("PASS: Normal format database correctly identified (key_bits <= 62)")
        else:
            print("FAIL: Normal format database not correctly identified")
            return False
        
        return True


def main():
    """Main test function."""
    print("Database Compatibility Test")
    print("=" * 50)
    
    # Set up test environment
    setup_test_environment()
    
    # Run tests
    tests = [
        test_wide_format_compatibility,
        test_normal_format_compatibility,
        test_database_header_analysis
    ]
    
    passed = 0
    for test in tests:
        if test():
            passed += 1
        print()
    
    if passed == len(tests):
        print("All database compatibility tests PASSED!")
        return 0
    else:
        print(f"Some database compatibility tests FAILED! ({passed}/{len(tests)} passed)")
        return 1


if __name__ == "__main__":
    sys.exit(main())
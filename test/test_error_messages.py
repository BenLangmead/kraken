#!/usr/bin/env python3

"""
Author: Ben Langmead
Date: 2025-09-04

Test script to verify error messages for database compatibility issues.
"""

import sys
import tempfile
from pathlib import Path
from test_utils import (
    setup_test_environment,
    create_wide_format_database,
    create_normal_format_database,
    run_test_kraken,
    print_database_info
)


def test_wide_format_error_messages():
    """Test error messages when loading wide format databases."""
    print("Testing wide format database error messages...")
    
    # Create test sequences for k=32
    sequences = ['ACGTACGTACGTACGTACGTACGTACGTACGT', 'TGCATGCATGCATGCATGCATGCATGCATGCA']
    wide_db = create_wide_format_database(sequences, 32)
    
    # Test that the database loads successfully with current Kraken
    result = run_test_kraken()
    
    if result.returncode == 0:
        print("PASS: Wide format database loads successfully with current Kraken")
        return True
    else:
        print(f"FAIL: Wide format database failed to load: {result.stderr}")
        return False


def test_database_header_inspection():
    """Test that we can inspect database headers to understand the format."""
    print("Testing database header inspection...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Create wide format database
        wide_sequences = ['ACGTACGTACGTACGTACGTACGTACGTACGT', 'TGCATGCATGCATGCATGCATGCATGCATGCA']
        wide_db = create_wide_format_database(wide_sequences, 32, temp_dir)
        
        # Create normal format database
        normal_sequences = ['ACGTACGTACGTACGT', 'TGCATGCATGCATGCA']
        normal_db = create_normal_format_database(normal_sequences, 8, temp_dir)
        
        # Inspect headers
        print("\nDatabase Header Analysis:")
        print("-" * 40)
        
        print_database_info("Wide Format", wide_db)
        print()
        print_database_info("Normal Format", normal_db)
        print()
        
        return True


def main():
    """Main test function."""
    print("Database Error Message Test")
    print("=" * 50)
    
    # Set up test environment
    setup_test_environment()
    
    # Run tests
    tests = [
        test_wide_format_error_messages,
        test_database_header_inspection
    ]
    
    passed = 0
    for test in tests:
        if test():
            passed += 1
        print()
    
    if passed == len(tests):
        print("All error message tests PASSED!")
        return 0
    else:
        print(f"Some error message tests FAILED! ({passed}/{len(tests)} passed)")
        return 1


if __name__ == "__main__":
    sys.exit(main())
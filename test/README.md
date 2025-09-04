# Kraken Test Suite

This directory contains tests for the Kraken database building system.

## Running Tests

From the project root directory:

```bash
./test/run_tests.sh
```

Or run the Python driver directly:

```bash
python3 test/run_tests.py
```

## Test Components

- **`run_tests.py`**: Python driver script that discovers and runs all test modules
- **`run_tests.sh`**: Shell wrapper script for running tests from the project root

- **`test_kmc_conversion.py`**: Tests the KMC3 to Jellyfish format conversion tool
- **`test_integration.py`**: Integration test comparing Jellyfish vs KMC3 database building workflows

## Requirements

- Python 3
- Jellyfish v1.x
- KMC3
- Kraken build tools (compiled from source)




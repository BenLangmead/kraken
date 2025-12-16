#!/usr/bin/env python3

"""
Author: Ben Langmead
Date: 2025-09-04

Shared utilities for Kraken test scripts.
"""

import os
import sys
import subprocess
import tempfile
import shutil
from pathlib import Path


def run_command(cmd, check=True, capture_output=True):
    """Run a command and return the result."""
    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=capture_output, text=True)
    if check and result.returncode != 0:
        print(f"Command failed with return code {result.returncode}")
        if result.stdout:
            print(f"stdout: {result.stdout}")
        if result.stderr:
            print(f"stderr: {result.stderr}")
        sys.exit(1)
    return result


def create_test_fasta(filename, sequences):
    """Create a test FASTA file with given sequences."""
    with open(filename, 'w') as f:
        for i, seq in enumerate(sequences):
            f.write(f">test_{i+1}\n{seq}\n")


def create_test_multiline_fasta(filename, sequence):
    """Create a test FASTA file with given sequences."""
    with open(filename, 'w') as f:
        f.write(f">test_0\n")
        for j in range(0, len(sequence), 80):
            f.write(f"{sequence[j:j+80]}\n")

def check_required_tools(tools):
    """Check that all required tools are available in PATH."""
    for tool in tools:
        if not shutil.which(tool):
            raise RuntimeError(f"ERROR: Required tool '{tool}' not found in PATH")


def check_required_files(files):
    """Check that all required files exist."""
    for file_path in files:
        if not os.path.exists(file_path):
            raise RuntimeError(f"ERROR: Required file '{file_path}' not found")


def create_wide_format_database(sequences, k, multiline_1sequence, temp_dir=None):
    """Create a wide format database (k > 31) for testing."""
    if temp_dir is None:
        temp_dir = tempfile.mkdtemp()
    
    temp_path = Path(temp_dir)
    
    # Create test FASTA file
    fasta_file = temp_path / "test.fa"
    if multiline_1sequence:
        create_test_multiline_fasta(fasta_file, sequences[0])
    else:
        create_test_fasta(fasta_file, sequences)
    
    # Build KMC database
    kmc_prefix = temp_path / "kmc_db"
    kmc_cmd = f"kmc -k{k} -ci1 -fm {fasta_file} {kmc_prefix} {temp_path}"
    run_command(kmc_cmd)
    
    # Convert to wide format
    wide_db = temp_path / "wide.jdb"
    # Check if we're running from test directory or project root
    if os.path.exists('src/kmc_to_jellyfish'):
        # Running from project root
        convert_cmd = f"src/kmc_to_jellyfish -k {k} {kmc_prefix} {wide_db}"
    else:
        # Running from test directory
        convert_cmd = f"../src/kmc_to_jellyfish -k {k} {kmc_prefix} {wide_db}"
    run_command(convert_cmd)
    
    return wide_db


def create_normal_format_database(sequences, k, temp_dir=None):
    """Create a normal format database (k ≤ 31) for testing."""
    if temp_dir is None:
        temp_dir = tempfile.mkdtemp()
    
    temp_path = Path(temp_dir)
    
    # Create test FASTA file
    fasta_file = temp_path / "test.fa"
    create_test_fasta(fasta_file, sequences)
    
    # Build Jellyfish database
    jellyfish_prefix = temp_path / "jellyfish_db"
    jellyfish_cmd = f"jellyfish count -m {k} -s 1000 -C {fasta_file} -o {jellyfish_prefix}"
    run_command(jellyfish_cmd)
    
    # Handle multiple files
    jellyfish_files = list(temp_path.glob("jellyfish_db_*"))
    if len(jellyfish_files) > 1:
        jellyfish_merge_cmd = f"jellyfish merge -o {jellyfish_prefix}.jdb {' '.join(str(f) for f in jellyfish_files)}"
        run_command(jellyfish_merge_cmd)
    else:
        shutil.move(jellyfish_files[0], f"{jellyfish_prefix}.jdb")
    
    return f"{jellyfish_prefix}.jdb"


def analyze_database_header(db_path):
    """Analyze a database header and return key information."""
    with open(db_path, 'rb') as f:
        header = f.read(64)
    
    magic = header[:8]
    key_bits = int.from_bytes(header[8:16], byteorder='little')
    val_len = int.from_bytes(header[16:24], byteorder='little')
    key_len = int.from_bytes(header[24:32], byteorder='little')
    key_ct = int.from_bytes(header[48:56], byteorder='little')
    
    k = key_bits // 2
    format_type = 'Wide (k > 31)' if k > 31 else 'Normal (k ≤ 31)'
    
    return {
        'magic': magic,
        'key_bits': key_bits,
        'val_len': val_len,
        'key_len': key_len,
        'key_ct': key_ct,
        'k': k,
        'format_type': format_type
    }


def print_database_info(db_name, db_path):
    """Print formatted database information."""
    info = analyze_database_header(db_path)
    print(f"{db_name}:")
    print(f"  Magic: {info['magic']}")
    print(f"  Key bits: {info['key_bits']}")
    print(f"  K-mer length: {info['k']}")
    print(f"  Value length: {info['val_len']} bytes")
    print(f"  Key length: {info['key_len']} bytes")
    print(f"  K-mer count: {info['key_ct']}")
    print(f"  Format: {info['format_type']}")


def run_test_kraken():
    """Run the test_kraken program and return the result."""
    # Check if we're running from test directory or project root
    if os.path.exists('src/test_kraken'):
        # Running from project root
        test_cmd = "src/test_kraken"
    else:
        # Running from test directory
        test_cmd = "../src/test_kraken"
    return run_command(test_cmd, check=False, capture_output=True)


def setup_test_environment():
    """Set up the test environment by checking required tools and files."""
    required_tools = ['kmc', 'kmc_tools', 'jellyfish']
    
    # Check if we're running from test directory or project root
    if os.path.exists('src/kmc_to_jellyfish'):
        # Running from project root
        required_files = ['src/kmc_to_jellyfish', 'src/test_kraken']
    else:
        # Running from test directory
        required_files = ['../src/kmc_to_jellyfish', '../src/test_kraken']
    
    check_required_tools(required_tools)
    check_required_files(required_files)

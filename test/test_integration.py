#!/usr/bin/env python3

"""
Integration test for Kraken database building with Jellyfish vs KMC+conversion.

This script tests the k-mer counting step (step 1) to validate that KMC+conversion 
produces the same results as Jellyfish. It focuses on the core functionality needed
to ensure the KMC integration works correctly.
"""

import os
import sys
import subprocess
import tempfile
import time
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional

def run_command(cmd: str, check: bool = True) -> subprocess.CompletedProcess:
    """Run a shell command and return the result."""
    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if check and result.returncode != 0:
        print(f"Command failed with exit code {result.returncode}")
        print(f"STDOUT: {result.stdout}")
        print(f"STDERR: {result.stderr}")
        sys.exit(1)
    return result

def create_test_dataset(temp_dir: Path) -> Path:
    """Create a small test dataset with a few bacterial genomes."""
    library_dir = temp_dir / "library"
    library_dir.mkdir()
    
    # Create a few small bacterial genomes for testing
    genomes = [
        ("E_coli", "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"),
        ("S_aureus", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
        ("P_aeruginosa", "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA")
    ]
    
    for name, sequence in genomes:
        genome_file = library_dir / f"{name}.fna"
        with open(genome_file, 'w') as f:
            f.write(f">{name}\n{sequence}\n")
    
    print(f"Created test dataset with {len(genomes)} genomes")
    return library_dir

def build_jellyfish_database(temp_dir: Path, library_dir: Path, 
                           kmer_len: int = 31, threads: int = 1) -> Tuple[Path, float]:
    """Build a Jellyfish database using the original method."""
    print("\n=== Building Jellyfish database ===")
    
    # Set up environment variables
    env = os.environ.copy()
    env.update({
        'KRAKEN_DB_NAME': str(temp_dir / "jellyfish_db"),
        'KRAKEN_KMER_LEN': str(kmer_len),
        'KRAKEN_THREAD_CT': str(threads),
        'KRAKEN_HASH_SIZE': '1000',  # Small hash size for testing
        'KRAKEN_MINIMIZER_LEN': '15',
        'KRAKEN_MAX_DB_SIZE': '',  # No size limit for testing
        'KRAKEN_REBUILD_DATABASE': '1',
        'KRAKEN_WORK_ON_DISK': '1'  # Minimize RAM usage
    })
    
    # Create database directory
    db_dir = Path(env['KRAKEN_DB_NAME'])
    db_dir.mkdir(exist_ok=True)
    
    # Create library directory
    (db_dir / "library").mkdir(exist_ok=True)
    
    # Copy test genomes
    shutil.copytree(library_dir, db_dir / "library", dirs_exist_ok=True)
    
    start_time = time.time()
    
    # Run only the k-mer counting step (step 1)
    print("Running Jellyfish k-mer counting...")
    
    # Check for Jellyfish
    run_command("bash ../scripts/check_for_jellyfish.sh", check=True)
    
    # Count k-mers with Jellyfish
    jellyfish_cmd = f"find {db_dir}/library/ -name '*.fna' -print0 | xargs -0 cat | jellyfish count -m {kmer_len} -s 1000 -C -t {threads} -o {db_dir}/database /dev/fd/0"
    run_command(jellyfish_cmd)
    
    # Merge only if necessary
    if os.path.exists(f"{db_dir}/database_1"):
        merge_cmd = f"jellyfish merge -o {db_dir}/database.jdb {db_dir}/database_*"
        run_command(merge_cmd)
    else:
        shutil.move(f"{db_dir}/database_0", f"{db_dir}/database.jdb")
    
    build_time = time.time() - start_time
    
    print(f"Jellyfish database built in {build_time:.2f} seconds")
    return db_dir, build_time

def build_kmc_database(temp_dir: Path, library_dir: Path,
                      kmer_len: int = 31, threads: int = 1) -> Tuple[Path, float]:
    """Build a database using KMC+conversion."""
    print("\n=== Building KMC database ===")
    
    # Set up environment variables
    env = os.environ.copy()
    env.update({
        'KRAKEN_DB_NAME': str(temp_dir / "kmc_db"),
        'KRAKEN_KMER_LEN': str(kmer_len),
        'KRAKEN_THREAD_CT': str(threads),
        'KRAKEN_HASH_SIZE': '1000',  # Small hash size for testing
        'KRAKEN_MINIMIZER_LEN': '15',
        'KRAKEN_MAX_DB_SIZE': '',  # No size limit for testing
        'KRAKEN_REBUILD_DATABASE': '1',
        'KRAKEN_WORK_ON_DISK': '1'  # Minimize RAM usage
    })
    
    # Add src directory to PATH for Kraken tools
    src_path = os.path.abspath('../src')
    env['PATH'] = f"{src_path}:{env.get('PATH', '')}"
    
    # Create database directory
    db_dir = Path(env['KRAKEN_DB_NAME'])
    db_dir.mkdir(exist_ok=True)
    
    # Create library directory
    (db_dir / "library").mkdir(exist_ok=True)
    
    # Copy test genomes
    shutil.copytree(library_dir, db_dir / "library", dirs_exist_ok=True)
    
    start_time = time.time()
    
    # Run KMC k-mer counting
    print("Running KMC k-mer counting...")
    
    # Check for KMC tools (now integrated into check_for_jellyfish.sh)
    env['KRAKEN_USE_KMC'] = '1'
    run_command("bash ../scripts/check_for_jellyfish.sh", check=True)
    
    # Create temporary directory for KMC
    kmc_temp_dir = temp_dir / "kmc_temp"
    kmc_temp_dir.mkdir(exist_ok=True)
    
    # Create a combined FASTA file for KMC
    combined_fasta = db_dir / "combined.fna"
    with open(combined_fasta, 'w') as f:
        for genome_file in db_dir.glob("library/*.fna"):
            with open(genome_file, 'r') as genome:
                f.write(genome.read())
    
    # Count k-mers with KMC
    kmc_cmd = f"kmc -k{kmer_len} -fa {combined_fasta} {db_dir}/database_kmc {kmc_temp_dir}"
    run_command(kmc_cmd)
    
    # Convert KMC output to Jellyfish format
    convert_cmd = f"{src_path}/kmc_to_jellyfish -k {kmer_len} {db_dir}/database_kmc {db_dir}/database.jdb"
    run_command(convert_cmd)
    
    # Clean up KMC temporary files
    shutil.rmtree(kmc_temp_dir)
    for ext in ['.kmc_pre', '.kmc_suf']:
        if os.path.exists(f"{db_dir}/database_kmc{ext}"):
            os.remove(f"{db_dir}/database_kmc{ext}")
    
    build_time = time.time() - start_time
    
    print(f"KMC database built in {build_time:.2f} seconds")
    return db_dir, build_time

def compare_databases(jellyfish_db_dir: Path, kmc_db_dir: Path) -> Dict[str, bool]:
    """Compare the databases built with different tools."""
    print("\n=== Comparing Databases ===")
    
    results = {}
    
    # Compare database.jdb files
    jellyfish_db = jellyfish_db_dir / "database.jdb"
    kmc_db = kmc_db_dir / "database.jdb"
    
    if jellyfish_db.exists() and kmc_db.exists():
        # Compare file sizes
        jellyfish_size = jellyfish_db.stat().st_size
        kmc_size = kmc_db.stat().st_size
        results['size_match'] = jellyfish_size == kmc_size
        
        print(f"Database file sizes: Jellyfish={jellyfish_size}, KMC={kmc_size}")
        
        # Compare k-mer content using jellyfish dump
        try:
            jellyfish_dump = subprocess.run(f"jellyfish dump {jellyfish_db}", 
                                          shell=True, capture_output=True, text=True)
            kmc_dump = subprocess.run(f"jellyfish dump {kmc_db}", 
                                    shell=True, capture_output=True, text=True)
            
            if jellyfish_dump.returncode == 0 and kmc_dump.returncode == 0:
                results['content_match'] = jellyfish_dump.stdout == kmc_dump.stdout
                print(f"K-mer content: {'MATCH' if results['content_match'] else 'DIFFERENT'}")
                
                if not results['content_match']:
                    print("Jellyfish dump output:")
                    print(jellyfish_dump.stdout[:500] + "..." if len(jellyfish_dump.stdout) > 500 else jellyfish_dump.stdout)
                    print("KMC dump output:")
                    print(kmc_dump.stdout[:500] + "..." if len(kmc_dump.stdout) > 500 else kmc_dump.stdout)
            else:
                results['content_match'] = False
                print("Failed to dump databases for comparison")
        except Exception as e:
            results['content_match'] = False
            print(f"Error comparing k-mer content: {e}")
    else:
        results['size_match'] = False
        results['content_match'] = False
        print("One or both database files missing")
    
    return results

def main():
    """Main integration test function."""
    print("Kraken Database Building Integration Test")
    print("=" * 60)
    
    # Check required tools
    required_tools = ['jellyfish', 'kmc', 'kmc_tools']
    for tool in required_tools:
        if shutil.which(tool) is None:
            print(f"ERROR: Required tool '{tool}' not found in PATH")
            sys.exit(1)
    
    # Check if kmc_to_jellyfish tool exists
    if not os.path.exists("../src/kmc_to_jellyfish"):
        print("ERROR: kmc_to_jellyfish tool not found. Please build it first with 'make -C src'")
        sys.exit(1)
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Create test dataset
        library_dir = create_test_dataset(temp_path)
        
        # Build databases with both methods
        jellyfish_db_dir, jellyfish_time = build_jellyfish_database(
            temp_path, library_dir, kmer_len=31, threads=1)
        
        kmc_db_dir, kmc_time = build_kmc_database(
            temp_path, library_dir, kmer_len=31, threads=1)
        
        # Compare databases
        db_comparison = compare_databases(jellyfish_db_dir, kmc_db_dir)
        
        # Report results
        print("\n" + "=" * 60)
        print("INTEGRATION TEST RESULTS")
        print("=" * 60)
        
        print(f"Build times:")
        print(f"  Jellyfish: {jellyfish_time:.2f} seconds")
        print(f"  KMC: {kmc_time:.2f} seconds")
        print(f"  Speedup: {jellyfish_time/kmc_time:.2f}x")
        
        print(f"\nDatabase comparison:")
        print(f"  File size match: {'PASS' if db_comparison['size_match'] else 'FAIL'}")
        print(f"  K-mer content match: {'PASS' if db_comparison['content_match'] else 'FAIL'}")
        
        # Overall result
        all_passed = (db_comparison['size_match'] and 
                     db_comparison['content_match'])
        
        print(f"\nOverall result: {'PASS' if all_passed else 'FAIL'}")
        
        if all_passed:
            print("All integration tests PASSED!")
            sys.exit(0)
        else:
            print("Some integration tests FAILED")
            sys.exit(1)

if __name__ == "__main__":
    main()

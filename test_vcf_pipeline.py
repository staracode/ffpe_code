#!/usr/bin/env python3
"""
Unit tests for FFPE VCF processing pipeline.

This test suite validates the complete workflow:
1. VCF format conversion
2. Mutation counting and analysis
3. Output validation

Test data is located in test_data/ directory.
"""

import unittest
import os
import tempfile
import shutil
import subprocess
import sys
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'scripts'))

# Try to import required modules, with graceful handling of missing dependencies
try:
    from convert_vcf import convert_vcf
    CONVERT_VCF_AVAILABLE = True
except ImportError as e:
    print(f"Warning: convert_vcf module not available: {e}")
    CONVERT_VCF_AVAILABLE = False

try:
    from mutation import count_mutations
    MUTATION_AVAILABLE = True
except ImportError as e:
    print(f"Warning: mutation module not available (missing dependencies): {e}")
    print("Required dependencies: cyvcf2, pyfaidx, biopython")
    print("Install with: conda install -c bioconda cyvcf2 pyfaidx biopython")
    MUTATION_AVAILABLE = False


class TestVCFPipeline(unittest.TestCase):
    """Test suite for VCF processing pipeline."""
    
    def setUp(self):
        """Set up test environment."""
        self.test_data_dir = Path("test_data")
        self.temp_dir = tempfile.mkdtemp()
        self.genome_dir = Path("genome_files")
        
        # Test file paths
        self.input_vcf = self.test_data_dir / "study3_FFPE-only_mutations.vcf"
        self.expected_adj_vcf = self.test_data_dir / "study3_FFPE-only_mutations_adj.vcf"
        self.expected_counts = self.test_data_dir / "study3_FFPE-only_mutations_count.csv"
        
        # Temporary output files
        self.temp_adj_vcf = Path(self.temp_dir) / "test_adj.vcf"
        self.temp_counts = Path(self.temp_dir) / "test_counts.tsv"
        
        # Reference genome (assuming it exists)
        self.reference_fasta = Path(self.genome_dir)/"hg19.fa"
        
        # Verify test data exists
        self.assertTrue(self.input_vcf.exists(), f"Test input VCF not found: {self.input_vcf}")
        self.assertTrue(self.expected_adj_vcf.exists(), f"Expected adjusted VCF not found: {self.expected_adj_vcf}")
        self.assertTrue(self.expected_counts.exists(), f"Expected counts file not found: {self.expected_counts}")
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    @unittest.skipUnless(CONVERT_VCF_AVAILABLE, "convert_vcf module not available")
    def test_vcf_conversion(self):
        """Test VCF format conversion."""
        print(f"\nTesting VCF conversion: {self.input_vcf} -> {self.temp_adj_vcf}")
        
        # Convert VCF
        convert_vcf(str(self.input_vcf), str(self.temp_adj_vcf))
        
        # Verify output file was created
        self.assertTrue(self.temp_adj_vcf.exists(), "Converted VCF file was not created")
        
        # Read and validate converted VCF
        with open(self.temp_adj_vcf, 'r') as f:
            lines = f.readlines()
        
        # Check header
        self.assertTrue(lines[0].startswith("##fileformat=VCFv4.2"), "Missing VCF header")
        self.assertTrue(lines[1].startswith("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), "Missing column header")
        
        # Check data lines
        data_lines = [line for line in lines if not line.startswith("#")]
        self.assertGreater(len(data_lines), 0, "No data lines found in converted VCF")
        
        # Validate format of first few data lines
        for line in data_lines[:5]:
            parts = line.strip().split("\t")
            self.assertEqual(len(parts), 8, f"Invalid number of columns: {len(parts)}")
            self.assertIsInstance(parts[1], str, "Position should be string")
            self.assertTrue(parts[1].isdigit(), "Position should be numeric")
        
        print(f"✓ VCF conversion successful: {len(data_lines)} variants processed")
    
    @unittest.skipUnless(MUTATION_AVAILABLE, "mutation module not available (missing dependencies)")
    def test_mutation_counting(self):
        """Test mutation counting from VCF."""
        if not os.path.exists(self.reference_fasta):
            self.skipTest(f"Reference genome not found: {self.reference_fasta}")
        
        print(f"\nTesting mutation counting: {self.expected_adj_vcf} -> {self.temp_counts}")
        
        # Count mutations
        count_mutations(str(self.expected_adj_vcf), self.reference_fasta, str(self.temp_counts))
        
        # Verify output file was created
        self.assertTrue(self.temp_counts.exists(), "Mutation counts file was not created")
        
        # Read and validate counts
        with open(self.temp_counts, 'r') as f:
            lines = f.readlines()
        
        # Check header
        self.assertTrue(lines[0].startswith("MutationType\tCount"), "Missing counts header")
        
        # Check data
        data_lines = [line for line in lines[1:] if line.strip()]
        self.assertGreater(len(data_lines), 0, "No mutation counts found")
        
        # Validate format and content
        total_mutations = 0
        for line in data_lines:
            parts = line.strip().split("\t")
            self.assertEqual(len(parts), 2, f"Invalid mutation count format: {line}")
            
            mutation_type, count = parts
            self.assertTrue(count.isdigit(), f"Count should be numeric: {count}")
            count_val = int(count)
            total_mutations += count_val
            
            # Validate mutation type format (e.g., "C>A@A_A")
            self.assertIn("@", mutation_type, f"Invalid mutation type format: {mutation_type}")
        
        print(f"✓ Mutation counting successful: {len(data_lines)} mutation types, {total_mutations} total mutations")
    
    def test_complete_workflow(self):
        """Test the complete workflow: VCF conversion -> mutation counting."""
        if not os.path.exists(self.reference_fasta):
            self.skipTest(f"Reference genome not found: {self.reference_fasta}")
        
        print(f"\nTesting complete workflow...")
        
        # Step 1: Convert VCF
        convert_vcf(str(self.input_vcf), str(self.temp_adj_vcf))
        self.assertTrue(self.temp_adj_vcf.exists())
        
        # Step 2: Count mutations
        count_mutations(str(self.temp_adj_vcf), self.reference_fasta, str(self.temp_counts))
        self.assertTrue(self.temp_counts.exists())
        
        # Step 3: Validate output against expected results
        self._compare_with_expected_results()
        
        print("✓ Complete workflow test passed")
    
    def test_command_line_interface(self):
        """Test command line interface of the scripts."""
        print(f"\nTesting command line interfaces...")
        
        # Test convert_vcf.py CLI
        try:
            result = subprocess.run([
                sys.executable, "scripts/convert_vcf.py", 
                str(self.input_vcf), str(self.temp_adj_vcf)
            ], capture_output=True, text=True, check=True)
            self.assertTrue(self.temp_adj_vcf.exists())
            print("✓ convert_vcf.py CLI test passed")
        except subprocess.CalledProcessError as e:
            self.fail(f"convert_vcf.py CLI failed: {e.stderr}")
        
        # Test mutation.py CLI (if reference exists)
        if os.path.exists(self.reference_fasta):
            try:
                result = subprocess.run([
                    sys.executable, "scripts/mutation.py",
                    str(self.temp_adj_vcf), self.reference_fasta, str(self.temp_counts)
                ], capture_output=True, text=True, check=True)
                self.assertTrue(self.temp_counts.exists())
                print("✓ mutation.py CLI test passed")
            except subprocess.CalledProcessError as e:
                self.fail(f"mutation.py CLI failed: {e.stderr}")
        else:
            print("⚠ Skipping mutation counter CLI test (reference genome not found)")
    
    def _compare_with_expected_results(self):
        """Compare generated results with expected results."""
        if not self.expected_counts.exists():
            self.skipTest("Expected results file not found")
        
        # Read expected results
        expected_counts = {}
        with open(self.expected_counts, 'r') as f:
            next(f)  # Skip header
            for line in f:
                mutation_type, count = line.strip().split("\t")
                expected_counts[mutation_type] = int(count)
        
        # Read generated results
        generated_counts = {}
        with open(self.temp_counts, 'r') as f:
            next(f)  # Skip header
            for line in f:
                mutation_type, count = line.strip().split("\t")
                generated_counts[mutation_type] = int(count)
        
        # Compare results
        self.assertEqual(set(expected_counts.keys()), set(generated_counts.keys()), 
                        "Mutation types don't match")
        
        differences = []
        for mutation_type in expected_counts:
            expected = expected_counts[mutation_type]
            generated = generated_counts[mutation_type]
            if expected != generated:
                differences.append(f"{mutation_type}: expected {expected}, got {generated}")
        
        if differences:
            self.fail(f"Count differences found:\n" + "\n".join(differences))
        
        print("✓ Results match expected output exactly")


class TestVCFValidation(unittest.TestCase):
    """Additional validation tests for VCF processing."""
    
    def test_vcf_format_validation(self):
        """Test VCF format validation."""
        from convert_vcf import convert_vcf
        
        # Create a minimal test VCF
        test_vcf_content = """#Chr\tpos\tanno\tref\talt\tRD_ref\tRD_alt
chr1\t1000\t.\tA\tT\t100\t10
chr2\t2000\t.\tC\tG\t150\t5"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(test_vcf_content)
            input_file = f.name
        
        output_file = input_file.replace('.vcf', '_adj.vcf')
        
        try:
            convert_vcf(input_file, output_file)
            
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            # Validate format
            self.assertTrue(lines[0].startswith("##fileformat=VCFv4.2"))
            self.assertTrue(lines[1].startswith("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"))
            
            # Check data lines
            data_lines = [line for line in lines if not line.startswith("#")]
            self.assertEqual(len(data_lines), 2)
            
            for line in data_lines:
                parts = line.strip().split("\t")
                self.assertEqual(len(parts), 8)
                self.assertEqual(parts[2], ".")  # ID should be "."
                self.assertEqual(parts[5], ".")  # QUAL should be "."
                self.assertEqual(parts[6], ".")  # FILTER should be "."
                self.assertEqual(parts[7], ".")  # INFO should be "."
        
        finally:
            os.unlink(input_file)
            if os.path.exists(output_file):
                os.unlink(output_file)


def run_tests():
    """Run all tests with verbose output."""
    print("🧪 Running FFPE VCF Pipeline Tests")
    print("=" * 50)
    
    # Print dependency status
    print(f"✓ convert_vcf module: {'Available' if CONVERT_VCF_AVAILABLE else 'Not available'}")
    print(f"✓ mutation module: {'Available' if MUTATION_AVAILABLE else 'Not available (missing dependencies)'}")
    
    if not MUTATION_AVAILABLE:
        print("\n📦 To install missing dependencies:")
        print("conda install -c bioconda cyvcf2 pyfaidx biopython")
    
    print("=" * 50)
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestVCFPipeline))
    suite.addTests(loader.loadTestsFromTestCase(TestVCFValidation))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Print summary
    print("\n" + "=" * 50)
    if result.wasSuccessful():
        print("🎉 All tests passed!")
    else:
        print(f"❌ {len(result.failures)} failures, {len(result.errors)} errors")
    
    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1) 
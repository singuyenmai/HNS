import unittest
from pathlib import Path
import tempfile
import shutil
import sys

# Ensure src is in sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from modules.reporter import MergedReporter

class TestMergedReporter(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
        
    def test_report_success(self):
        tsv_hits = [
            {
                "Genome_ID": "s1",
                "Contig_ID": "c1",
                "Start": "10",
                "Stop": "200",
                "Strand": "+",
                "Bit_Score": "30.0",
                "Sequence_Length": "6"
            }
        ]
        
        out_tsv = Path(self.temp_dir) / "all.tsv"
        
        reporter = MergedReporter()
        success = reporter.run(tsv_hits, out_tsv)
        
        self.assertTrue(success)
        self.assertTrue(out_tsv.is_file())
        
        # Verify TSV contents
        tsv_lines = out_tsv.read_text().splitlines()
        self.assertEqual(len(tsv_lines), 2)
        self.assertEqual(tsv_lines[0], "Genome_ID\tContig_ID\tStart\tStop\tStrand\tBit_Score\tSequence_Length\tPredicted_protein\tLocation\tStpA_cov\tStpA_identity\tHNS_cov\tHNS_identity")
        self.assertEqual(tsv_lines[1], "s1\tc1\t10\t200\t+\t30.0\t6\tUnknown\tUnknown\t0\t0\t0\t0")

if __name__ == "__main__":
    unittest.main()

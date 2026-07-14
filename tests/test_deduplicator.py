import unittest
from unittest.mock import patch, MagicMock
from pathlib import Path
import tempfile
import shutil
import sys

# Ensure src is in sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from modules.deduplicator import CDHitDeduplicator

class TestCDHitDeduplicator(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
        
    @patch("subprocess.run")
    def test_deduplicate_success(self, mock_run):
        # Mock successful subprocess execution
        mock_result = MagicMock()
        mock_result.stdout = "CD-HIT run finished"
        mock_run.return_value = mock_result
        
        in_fasta = Path(self.temp_dir) / "in.faa"
        out_fasta = Path(self.temp_dir) / "out.faa"
        in_fasta.write_text(">seq\nMSEALK\n")
        
        dedup = CDHitDeduplicator()
        success = dedup.run(
            input_fasta=in_fasta,
            output_unique_fasta=out_fasta,
            threads=8,
            mem_mb=8000
        )
        
        self.assertTrue(success)
        mock_run.assert_called_once()
        args = mock_run.call_args[0][0]
        self.assertIn("apptainer", args)
        self.assertIn("cd-hit", args)
        self.assertIn("-T", args)
        self.assertIn("8", args)
        self.assertIn("-M", args)
        self.assertIn("8000", args)

if __name__ == "__main__":
    unittest.main()

import unittest
from unittest.mock import patch, MagicMock
from pathlib import Path
import tempfile
import gzip
import shutil
import sys

# Ensure src is in sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from modules.downloader import HMMDownloader

class TestHMMDownloader(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.output_path = Path(self.temp_dir) / "PF00816.hmm"
        
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
        
    @patch("urllib.request.urlopen")
    def test_download_success(self, mock_urlopen):
        # Prepare mock response data (gzipped mock HMM)
        raw_content = b"mock hmm content"
        compressed_content = gzip.compress(raw_content)
        
        mock_response = MagicMock()
        mock_response.read.return_value = compressed_content
        mock_urlopen.return_value.__enter__.return_value = mock_response
        
        downloader = HMMDownloader()
        res_path = downloader.run(str(self.output_path))
        
        self.assertEqual(res_path, self.output_path)
        self.assertTrue(self.output_path.is_file())
        with open(self.output_path, "rb") as f:
            self.assertEqual(f.read(), raw_content)
            
    def test_existing_file(self):
        # Create a fake file
        self.output_path.write_text("existing hmm")
        downloader = HMMDownloader()
        res_path = downloader.run(str(self.output_path))
        
        self.assertEqual(res_path, self.output_path)
        with open(self.output_path, "r") as f:
            self.assertEqual(f.read(), "existing hmm")

if __name__ == "__main__":
    unittest.main()

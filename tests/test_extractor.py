import unittest
from unittest.mock import patch, MagicMock
from pathlib import Path
import tempfile
import shutil
import sys

# Ensure src is in sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from modules.extractor import process_genome_worker, HNSScanner

class TestHNSScanner(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
        
    @patch("pyhmmer.plan7.HMMFile")
    @patch("pyrodigal.GeneFinder")
    @patch("pyhmmer.hmmsearch")
    def test_process_genome_worker(self, mock_hmmsearch, mock_genefinder, mock_hmmfile):
        # Setup mock HMM
        import pyhmmer
        mock_hmm = MagicMock()
        mock_hmm.alphabet = pyhmmer.easel.Alphabet.amino()
        mock_hmm_instance = mock_hmmfile.return_value.__enter__.return_value
        mock_hmm_instance.read.return_value = mock_hmm
        
        # Setup mock pyrodigal gene finder
        mock_gf_instance = mock_genefinder.return_value
        mock_gene = MagicMock()
        mock_gene.begin = 10
        mock_gene.end = 200
        mock_gene.strand = 1
        mock_gene.translate.return_value = "MSEALK"
        mock_gf_instance.find_genes.return_value = [mock_gene]
        
        # Setup mock pyhmmer search hits
        mock_hit = MagicMock()
        mock_hit.score = 30.0
        mock_hit.name = b"sample1||contig1||10||200||+||1"
        
        # hmmsearch yields a list/iterator of hits
        mock_hmmsearch.return_value = [[mock_hit]]
        
        # Create a tiny mock FASTA file
        fasta_file = Path(self.temp_dir) / "genome.fasta"
        fasta_file.write_text(">contig1\nATGCGTACGT\n")
        
        # Call process_genome_worker
        err, fasta_hits, tsv_hits, is_cached = process_genome_worker(
            sample_name="sample1",
            fasta_path_str=str(fasta_file),
            hmm_path_str="fake_hmm.hmm",
            bit_score_threshold=25.0,
            out_dir=self.temp_dir,
            hns_ref="MSEALK",
            stpa_ref="MSEALK"
        )
        
        self.assertIsNone(err)
        self.assertEqual(len(fasta_hits), 1)
        self.assertEqual(fasta_hits[0]["id"], "sample1||contig1||10||200||+||1")
        self.assertEqual(fasta_hits[0]["seq"], "MSEALK")
        
        self.assertEqual(len(tsv_hits), 1)
        self.assertEqual(tsv_hits[0]["Genome_ID"], "sample1")
        self.assertEqual(tsv_hits[0]["Contig_ID"], "contig1")
        self.assertEqual(tsv_hits[0]["Start"], "10")
        self.assertEqual(tsv_hits[0]["Stop"], "200")
        self.assertEqual(tsv_hits[0]["Strand"], "+")
        self.assertEqual(tsv_hits[0]["Bit_Score"], "30.0")
        self.assertEqual(tsv_hits[0]["Sequence_Length"], "6")
        
        # Check files are created
        self.assertTrue((Path(self.temp_dir) / "sample1.faa").is_file())
        self.assertTrue((Path(self.temp_dir) / "sample1.contigs.fasta").is_file())
        self.assertTrue((Path(self.temp_dir) / "sample1.done").is_file())
        self.assertTrue((Path(self.temp_dir) / "sample1_hits.json").is_file())
        
        # Test resume behavior by calling process_genome_worker again.
        # It should bypass mocks and directly load hits from the JSON cache.
        mock_genefinder.reset_mock()
        mock_hmmfile.reset_mock()
        
        err2, fasta_hits2, tsv_hits2, is_cached2 = process_genome_worker(
            sample_name="sample1",
            fasta_path_str=str(fasta_file),
            hmm_path_str="fake_hmm.hmm",
            bit_score_threshold=25.0,
            out_dir=self.temp_dir,
            hns_ref="MSEALK",
            stpa_ref="MSEALK"
        )
        
        self.assertIsNone(err2)
        self.assertEqual(len(fasta_hits2), 1)
        self.assertEqual(fasta_hits2[0]["id"], "sample1||contig1||10||200||+||1")
        # Ensure mocks were not called again during resume
        mock_genefinder.assert_not_called()
        mock_hmmfile.assert_not_called()

if __name__ == "__main__":
    unittest.main()

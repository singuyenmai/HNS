import unittest
from unittest.mock import patch, MagicMock
from pathlib import Path
import tempfile
import shutil
import sys
import pandas as pd

# Ensure src is in sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from modules.classifier import PlatonClassifier

class TestPlatonClassifier(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
        
    @patch("subprocess.run")
    def test_platon_classifier_run(self, mock_run):
        # Create a mock tsv_hits
        tsv_hits = [
            {"Genome_ID": "sample1", "Contig_ID": "contig1", "Start": "10"},
            {"Genome_ID": "sample1", "Contig_ID": "contig2", "Start": "20"},
            {"Genome_ID": "sample2", "Contig_ID": "contig3", "Start": "30"}
        ]
        unique_genomes = [("sample1", "fake_path1"), ("sample2", "fake_path2")]
        
        # We need to simulate the contigs.fasta file existing so Platon is called
        out_path = Path(self.temp_dir)
        (out_path / "sample1.contigs.fasta").touch()
        (out_path / "sample2.contigs.fasta").touch()
        
        # We need to create a side effect for subprocess.run that writes out the fake platon TSV files
        def side_effect_run(cmd, *args, **kwargs):
            input_fasta = Path(cmd[-1])
            sample_name = input_fasta.name.split(".")[0]
            
            platon_dir = out_path / f"{sample_name}_platon"
            platon_dir.mkdir(exist_ok=True)
            if sample_name == "sample1":
                with open(platon_dir / f"{sample_name}.contigs.plasmid.fasta", "w") as f:
                    f.write(">contig1\nATGC\n")
                with open(platon_dir / f"{sample_name}.contigs.chromosome.fasta", "w") as f:
                    f.write(">contig2\nATGC\n")
            else:
                with open(platon_dir / f"{sample_name}.contigs.plasmid.fasta", "w") as f:
                    f.write(">contig3\nATGC\n")
            
            return MagicMock(returncode=0)
            
        mock_run.side_effect = side_effect_run
        
        classifier = PlatonClassifier()
        updated_hits = classifier.run(tsv_hits, unique_genomes, self.temp_dir, cpus=1)
        
        # Assertions
        self.assertEqual(len(updated_hits), 3)
        self.assertEqual(updated_hits[0]["Location"], "Plasmid")
        self.assertEqual(updated_hits[1]["Location"], "Chromosome")
        self.assertEqual(updated_hits[2]["Location"], "Plasmid")
        
        self.assertEqual(mock_run.call_count, 2)

if __name__ == "__main__":
    unittest.main()

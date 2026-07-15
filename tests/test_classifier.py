import unittest
from unittest.mock import patch, MagicMock
from pathlib import Path
import tempfile
import shutil
import sys
import pandas as pd

# Ensure src is in sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from modules.classifier import LocationClassifier

class TestLocationClassifier(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
        
    @patch("subprocess.run")
    def test_platon_engine(self, mock_run):
        tsv_hits = [
            {"Genome_ID": "sample1", "Contig_ID": "contig1", "Start": "10"},
            {"Genome_ID": "sample1", "Contig_ID": "contig2", "Start": "20"},
            {"Genome_ID": "sample2", "Contig_ID": "contig3", "Start": "30"}
        ]
        unique_genomes = [("sample1", "fake_path1"), ("sample2", "fake_path2")]
        
        out_path = Path(self.temp_dir)
        (out_path / "sample1.contigs.fasta").touch()
        (out_path / "sample2.contigs.fasta").touch()
        
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
            
            # Create dummy tsv file to simulate done checkpoint
            (platon_dir / f"{sample_name}.contigs.tsv").touch()
            return MagicMock(returncode=0)
            
        mock_run.side_effect = side_effect_run
        
        classifier = LocationClassifier()
        updated_hits = classifier.run(tsv_hits, unique_genomes, self.temp_dir, cpus=1, classifier_engine="platon")
        
        self.assertEqual(len(updated_hits), 3)
        self.assertEqual(updated_hits[0]["Location"], "Plasmid")
        self.assertEqual(updated_hits[1]["Location"], "Chromosome")
        self.assertEqual(updated_hits[2]["Location"], "Plasmid")
        self.assertEqual(mock_run.call_count, 2)

    @patch("subprocess.run")
    def test_mlplasmids_engine(self, mock_run):
        tsv_hits = [
            {"Genome_ID": "sample1", "Contig_ID": "contig1", "Start": "10"},
            {"Genome_ID": "sample1", "Contig_ID": "contig2", "Start": "20"},
            {"Genome_ID": "sample2", "Contig_ID": "contig3", "Start": "30"}
        ]
        unique_genomes = [("sample1", "fake_path1"), ("sample2", "fake_path2")]
        
        out_path = Path(self.temp_dir)
        (out_path / "sample1.contigs.fasta").touch()
        (out_path / "sample2.contigs.fasta").touch()
        
        def side_effect_run(cmd, *args, **kwargs):
            r_str = cmd[-1]
            import re
            file_match = re.search(r"write\.table\(.*,\s*file\s*=\s*'([^']+)'", r_str)
            if file_match:
                ml_tsv_path = Path(file_match.group(1))
                ml_tsv_path.parent.mkdir(exist_ok=True, parents=True)
                
                if "sample1" in str(ml_tsv_path):
                    df = pd.DataFrame([
                        {"Contig_name": "contig1", "Prediction": "Plasmid"},
                        {"Contig_name": "contig2", "Prediction": "Chromosome"}
                    ])
                else:
                    df = pd.DataFrame([
                        {"Contig_name": "contig3", "Prediction": "Plasmid"}
                    ])
                df.to_csv(ml_tsv_path, sep="\t", index=False)
            
            return MagicMock(returncode=0)
            
        mock_run.side_effect = side_effect_run
        
        classifier = LocationClassifier()
        updated_hits = classifier.run(tsv_hits, unique_genomes, self.temp_dir, cpus=1, classifier_engine="mlplasmids")
        
        self.assertEqual(len(updated_hits), 3)
        self.assertEqual(updated_hits[0]["Location"], "Plasmid")
        self.assertEqual(updated_hits[1]["Location"], "Chromosome")
        self.assertEqual(updated_hits[2]["Location"], "Plasmid")
        self.assertEqual(mock_run.call_count, 2)

    @patch("subprocess.run")
    def test_genomad_engine(self, mock_run):
        tsv_hits = [
            {"Genome_ID": "sample1", "Contig_ID": "contig1", "Start": "10", "Stop": "15"},
            {"Genome_ID": "sample1", "Contig_ID": "contig2", "Start": "20", "Stop": "25"},
            {"Genome_ID": "sample1", "Contig_ID": "contig3", "Start": "30", "Stop": "35"},
            {"Genome_ID": "sample2", "Contig_ID": "contig4", "Start": "40", "Stop": "45"}
        ]
        unique_genomes = [("sample1", "fake_path1"), ("sample2", "fake_path2")]
        
        out_path = Path(self.temp_dir)
        (out_path / "sample1.contigs.fasta").touch()
        (out_path / "sample2.contigs.fasta").touch()
        
        def side_effect_run(cmd, *args, **kwargs):
            genomad_out_path = Path([x for x in cmd if "_genomad" in x][0])
            sample_name = genomad_out_path.name.replace("_genomad", "")
            
            summary_dir = genomad_out_path / f"{sample_name}.contigs_summary"
            summary_dir.mkdir(exist_ok=True, parents=True)
            
            if sample_name == "sample1":
                plasmid_df = pd.DataFrame([{"seq_name": "contig1|provirus_1"}])
                virus_df = pd.DataFrame([
                    {"seq_name": "contig2|provirus_15_25", "topology": "Provirus", "coordinates": "15-25"},
                    {"seq_name": "contig3|provirus_40_50", "topology": "Provirus", "coordinates": "40-50"}
                ])
            else:
                plasmid_df = pd.DataFrame([{"seq_name": "contig4"}])
                virus_df = pd.DataFrame(columns=["seq_name", "topology", "coordinates"])
                
            plasmid_df.to_csv(summary_dir / f"{sample_name}.contigs_plasmid_summary.tsv", sep="\t", index=False)
            virus_df.to_csv(summary_dir / f"{sample_name}.contigs_virus_summary.tsv", sep="\t", index=False)
            
            return MagicMock(returncode=0)
            
        mock_run.side_effect = side_effect_run
        
        classifier = LocationClassifier()
        updated_hits = classifier.run(tsv_hits, unique_genomes, self.temp_dir, cpus=1, classifier_engine="genomad")
        
        self.assertEqual(len(updated_hits), 4)
        self.assertEqual(updated_hits[0]["Location"], "Plasmid")
        self.assertEqual(updated_hits[1]["Location"], "Virus")
        self.assertEqual(updated_hits[2]["Location"], "Chromosome")
        self.assertEqual(updated_hits[3]["Location"], "Plasmid")
        self.assertEqual(mock_run.call_count, 2)

if __name__ == "__main__":
    unittest.main()

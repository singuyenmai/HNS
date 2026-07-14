import gzip
import sys
import urllib.request
from pathlib import Path
from .base import PipelineStep

class HMMDownloader(PipelineStep):
    """
    Downloads and decompresses the H-NS HMM profile from the InterPro Pfam database.
    """
    def run(self, output_path: str = "PF00816.hmm") -> Path:
        out_file = Path(output_path)
        if out_file.is_file():
            print(f"HMM profile {out_file} already exists.")
            return out_file

        url = "https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00816/?annotation=hmm"
        print(f"Downloading H-NS HMM profile from {url}...")
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req) as response:
                data = response.read()
            if data.startswith(b'\x1f\x8b'):
                print("Decompressing gzipped HMM profile...")
                data = gzip.decompress(data)
            with open(out_file, "wb") as f:
                f.write(data)
            print(f"Saved HMM profile to {out_file}")
            return out_file
        except Exception as e:
            print(f"Error downloading HMM profile: {e}", file=sys.stderr)
            raise

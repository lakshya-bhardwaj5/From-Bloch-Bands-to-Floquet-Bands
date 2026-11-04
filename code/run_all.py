"""Run from any working directory: python code/run_all.py [--output-run DIR]."""
import argparse
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
os.environ.setdefault("OMP_NUM_THREADS","1")
import json
import sys
import time
import hashlib
import platform
from pathlib import Path
import numpy as np
import scipy
import matplotlib
import static_bands, single_band, wavepackets, multiband, convergence

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output-run",type=Path,help="Optional fresh runner result directory; writes completion.json there.")
    args=parser.parse_args()
    start=time.time()
    for module in [static_bands,single_band,wavepackets,multiband,convergence]:
        print(f"Running {module.__name__}",flush=True)
        module.main()
    root=Path(__file__).resolve().parents[1]
    # A record cannot contain a hash of itself; exclude an older record left
    # by a previous run before constructing hashes of the scientific outputs.
    outputs=sorted(list((root/"figures").glob("*"))+list((root/"results").glob("*")))
    outputs=[p for p in outputs if p.name!="execution_record.json"]
    record={"status":"completed","elapsed_seconds":time.time()-start,
            "versions":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"matplotlib":matplotlib.__version__},
            "argv":sys.argv,"randomness":"none","outputs_sha256":{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs if p.is_file()}}
    (root/"results"/"execution_record.json").write_text(json.dumps(record,indent=2)+"\n")
    if args.output_run:
        args.output_run.mkdir(parents=True,exist_ok=True)
        (args.output_run/"completion.json").write_text(json.dumps(record,indent=2)+"\n")

if __name__=="__main__":main()

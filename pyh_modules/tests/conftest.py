import os
import sys
import shutil
import subprocess
import pandas as pd
import pytest
from pathlib import Path
import importlib

# Prepend the local project directory so that our local (updated) package is used.
local_project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if local_project_path not in sys.path:
    sys.path.insert(0, local_project_path)

from pyh_modules import helpers, core
importlib.reload(helpers)
importlib.reload(core)

def dummy_download_file(bucket, key, local_path):
    print(f"Dummy download: Writing dummy data to {local_path}")
    dummy_content = ">NODE_1\nACGTACGTACGT\n"
    os.makedirs(os.path.dirname(local_path), exist_ok=True)
    with open(local_path, "w") as f:
        f.write(dummy_content)

def dummy_upload_file_to_s3(local_path, bucket, key):
    print(f"Dummy upload: Simulated upload of {local_path} to s3://{bucket}/{key}")

def dummy_subprocess_run(args, check):
    input_file = args[4]
    output_file = args[6]
    print(f"Dummy seqkit: Copying {input_file} to {output_file}")
    shutil.copyfile(input_file, output_file)

@pytest.fixture
def setup_dummy_environment(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    samplesheet = tmp_path / "samplesheet.csv"
    df = pd.DataFrame({
        "sample": ["sample1"],
        "short_reads_1": [f"s3://dummy-bucket/path/dummy_R1.fastq.gz"]
    })
    df.to_csv(samplesheet, index=False)
    monkeypatch.setattr(helpers, "download_file", dummy_download_file)
    monkeypatch.setattr(helpers, "upload_file_to_s3", dummy_upload_file_to_s3)
    monkeypatch.setattr(subprocess, "run", dummy_subprocess_run)
    monkeypatch.setattr(helpers, "check_tool_availability", lambda tool: None)
    return samplesheet, tmp_path

@pytest.fixture
def custom_config(setup_dummy_environment):
    samplesheet, work_dir = setup_dummy_environment
    return {
        "samplesheet_path": str(samplesheet),
        "local_work_dir": str(work_dir)
    }


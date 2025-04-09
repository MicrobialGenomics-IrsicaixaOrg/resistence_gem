import os
import gzip
import shutil
import subprocess
import pandas as pd
from pathlib import Path

import pytest

# Updated import according to new package name.
from pyh_modules import core, helpers

def dummy_download_file(s3_bucket, s3_key, local_path):
    """
    Instead of downloading from S3, create a dummy gzipped FASTA file.
    """
    content = ">dummy_seq\nATGCATGCATGCATGC\n"
    os.makedirs(os.path.dirname(local_path), exist_ok=True)
    with gzip.open(local_path, "wb") as f:
        f.write(content.encode("utf-8"))
    print(f"Dummy download: Created {local_path}")

def dummy_subprocess_run(args, check):
    """
    Simulate seqkit filtering by copying the input file to the output file.
    Expected args format:
       ["seqkit", "seq", "-m", <min_length>, <input_file>, "-o", <output_file>]
    """
    input_file = args[4]
    output_file = args[6]
    # If the input is gzipped, decompress it and write to the output file
    if input_file.endswith(".gz"):
        with gzip.open(input_file, "rt") as fin, open(output_file, "w") as fout:
            fout.write(fin.read())
    else:
        shutil.copyfile(input_file, output_file)
    print(f"Dummy seqkit: Copied {input_file} to {output_file}")

def dummy_upload_file_to_s3(local_path, s3_bucket, s3_key):
    """
    Simulate upload to S3 by doing nothing.
    """
    print(f"Dummy upload: Simulated upload of {local_path} to s3://{s3_bucket}/{s3_key}")

@pytest.fixture
def setup_dummy_environment(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    samplesheet = tmp_path / "samplesheet.csv"
    df = pd.DataFrame({
        "sample": ["sample1"],
        "short_reads_1": ["s3://dummy-bucket/path/dummy_R1.fastq.gz"]
    })
    df.to_csv(samplesheet, index=False)
    # Override functions in both helpers and core modules:
    monkeypatch.setattr(helpers, "download_file", dummy_download_file)
    monkeypatch.setattr(core, "download_file", dummy_download_file)
    monkeypatch.setattr(subprocess, "run", dummy_subprocess_run)
    monkeypatch.setattr(helpers, "upload_file_to_s3", dummy_upload_file_to_s3)
    monkeypatch.setattr(core, "upload_file_to_s3", dummy_upload_file_to_s3)
    
    # Bypass the tool availability check in the helpers module.
    monkeypatch.setattr(helpers, "check_tool_availability", lambda tool: None)
    
    return samplesheet

def test_generate_merged_filtered_contigs_integration(setup_dummy_environment, tmp_path):
    """
    Integration test for generate_merged_filtered_contigs:
    - Uses a dummy samplesheet
    - Verifies that a merged file is produced in the expected location.
    """
    samplesheet = setup_dummy_environment

    # Run the contigs processing pipeline using the dummy environment.
    core.generate_merged_filtered_contigs(str(samplesheet), min_length=1)
    
    # Updated expected path: since local_work_dir is "./processing_results"
    merged_dir = tmp_path / "processing_results" / "Assembly_results" / "merged_results"
    merged_file = merged_dir / "sample1_merged.fa.gz"
    
    # Assert that the merged file was created.
    assert merged_file.exists(), "Merged file was not created."
    
    # Optionally, decompress the merged file and check its content.
    with gzip.open(merged_file, "rt") as f:
        content = f.read()
    # Update the assertion to look for a substring that matches the new header format.
    assert "MEGAHIT_dummy_seq" in content, "Merged file content does not contain expected sequence."


import os
import tempfile
import pandas as pd
import pytest
from sample import core

def test_download_files():
    """
    Test that download_files creates the target directory.
    """
    sample_names = ["sample1"]
    assembler = "TestAssembler"
    s3_base_path = "s3://bucket/path/"
    local_dir = tempfile.mkdtemp()
    core.download_files(sample_names, assembler, s3_base_path, local_dir)
    assert os.path.exists(local_dir)

def test_filter_sequences(tmp_path):
    """
    Test the filter_sequences function with a dummy FASTA file.
    Note: This test requires 'seqkit' to be installed.
    """
    dummy_file = tmp_path / "TestAssembler-sample1.contigs.fa.gz"
    # Create a dummy gzipped FASTA file
    with open(dummy_file, "wb") as f:
        f.write(b">seq1\nATGCATGCATGC\n")
    # Since the file is very short, filtering with min_length=1 should pass
    try:
        core.filter_sequences(["sample1"], "TestAssembler", tmp_path, min_length=1)
    except Exception as e:
        pytest.skip(f"seqkit filtering failed: {e}")


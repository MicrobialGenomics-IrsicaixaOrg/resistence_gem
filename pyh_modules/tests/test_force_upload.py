import os
import pytest
from pyh_modules import core, helpers

def dummy_upload_failure(local_path, s3_bucket, s3_key):
    """
    Dummy upload function that simulates failure for files containing 'fail_me'
    in their filename.
    """
    if "fail_me" in local_path:
        raise Exception("Simulated upload failure for force-upload testing.")
    else:
        print(f"Uploaded {local_path} successfully.")

@pytest.fixture
def force_upload_config(tmp_path):
    """
    Returns a configuration dictionary for testing both assemblies and bins with force_upload=True.
    A dummy samplesheet is created with a sample name that includes 'fail_me' so that our dummy upload
    function will simulate a failure.
    """
    samplesheet = tmp_path / "samplesheet.csv"
    samplesheet.write_text("sample,short_reads_1\nfail_me,s3://dummy-bucket/path/dummy_R1.fastq.gz\n")
    config = {
        "samplesheet_path": str(samplesheet),
        "local_work_dir": str(tmp_path / "results"),
        "force_upload": True  # Force re-upload even if merged files exist locally.
    }
    return config

def test_force_upload_assemblies(tmp_path, monkeypatch, force_upload_config):
    """
    Test that process_assemblies() forces a re-upload when force_upload is True.
    The dummy merged contigs file is created in the expected directory and the dummy upload function
    will always fail for this file. The test passes if no unhandled exception is raised.
    """
    # Patch the upload functions in both helpers and core to simulate failure.
    monkeypatch.setattr(helpers, "upload_file_to_s3", dummy_upload_failure)
    monkeypatch.setattr(core, "upload_file_to_s3", dummy_upload_failure)
    # Patch download_file so that no real downloads occur.
    monkeypatch.setattr(helpers, "download_file", lambda bucket, key, path: None)

    # Create a dummy merged contigs file in the expected assemblies location.
    merged_dir = os.path.join(force_upload_config["local_work_dir"],
                              "Assembly_results", "merged_results")
    os.makedirs(merged_dir, exist_ok=True)
    merged_file = os.path.join(merged_dir, "fail_me_merged.fa.gz")
    with open(merged_file, "w") as f:
        f.write(">dummy_seq\nATGCATGC\n")
    
    # Run the assemblies processing pipeline. It should attempt to re-upload the merged file.
    try:
        core.generate_merged_filtered_contigs(
            force_upload_config["samplesheet_path"],
            min_length=1,
            config=force_upload_config
        )
    except Exception as e:
        pytest.fail(f"process_assemblies() raised an unexpected exception: {e}")

def test_force_upload_bins(tmp_path, monkeypatch, force_upload_config):
    """
    Test that process_bins() forces a re-upload when force_upload is True.
    The dummy merged bins file is created in the expected location and the dummy upload function
    simulates failure for this file. The test passes if no unhandled exception is raised.
    """
    # Patch the upload functions to simulate failure.
    monkeypatch.setattr(helpers, "upload_file_to_s3", dummy_upload_failure)
    monkeypatch.setattr(core, "upload_file_to_s3", dummy_upload_failure)
    # Patch download_file so that no real downloads occur.
    monkeypatch.setattr(helpers, "download_file", lambda bucket, key, path: None)

    # Create a dummy merged bins file in the expected bins location.
    merged_bins_dir = os.path.join(force_upload_config["local_work_dir"],
                                   "Binning_results", "merged_bins")
    os.makedirs(merged_bins_dir, exist_ok=True)
    merged_file = os.path.join(merged_bins_dir, "fail_me_merged.fa.gz")
    with open(merged_file, "w") as f:
        f.write(">dummy_seq\nATGCATGC\n")
    
    # Run the bins processing pipeline. It should try to re-upload the merged bins file.
    try:
        core.generate_merged_filtered_bins(
            min_completeness=50,
            max_contamination=10,
            min_length=1,
            config=force_upload_config
        )
    except Exception as e:
        pytest.fail(f"process_bins() raised an unexpected exception: {e}")


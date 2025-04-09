import os
import pytest
from pyh_modules import core, helpers
import shutil

def failing_upload_file(local_path, s3_bucket, s3_key):
    """
    Dummy upload function that simulates an upload failure for files that contain 'fail_me'
    in their filename.
    """
    if "fail_me" in local_path:
        raise Exception("Simulated upload failure for testing force upload.")
    else:
        print(f"Uploaded {local_path} to s3://{s3_bucket}/{s3_key}")

@pytest.fixture
def force_upload_config(tmp_path):
    """
    Creates a dummy samplesheet and returns a configuration dictionary with force_upload=True.
    The samplesheet contains a sample name that includes 'fail_me' to trigger a simulated failure.
    """
    samplesheet = tmp_path / "samplesheet.csv"
    samplesheet.write_text("sample,short_reads_1\nfail_me,s3://dummy-bucket/path/dummy_R1.fastq.gz\n")
    config = {
        "samplesheet_path": str(samplesheet),
        "local_work_dir": str(tmp_path / "results"),
        "force_upload": True  # Enable forced re-upload
    }
    return config

def test_force_upload_behavior(tmp_path, monkeypatch, force_upload_config):
    """
    Tests that when force_upload=True, the process_bins pipeline attempts to upload the merged bin file,
    even if it already exists locally.
    
    The dummy upload function is set to fail when the filename contains 'fail_me', so we expect the error
    to be logged but not to crash the pipeline.
    """
    # Override the upload function to simulate failure for files containing "fail_me"
    monkeypatch.setattr(helpers, "upload_file_to_s3", failing_upload_file)
    monkeypatch.setattr(core, "upload_file_to_s3", failing_upload_file)
    
    # Override other external functions as needed so that they do not trigger network calls.
    monkeypatch.setattr(helpers, "download_file", lambda bucket, key, local: None)
    
    # Ensure the force_upload flag is enabled from the configuration.
    config = force_upload_config
    
    # Simulate that the merged bin file already exists.
    merged_bins_dir = os.path.join(config["local_work_dir"], "Binning_results", "merged_bins")
    os.makedirs(merged_bins_dir, exist_ok=True)
    merged_file_path = os.path.join(merged_bins_dir, "fail_me_merged.fa.gz")
    with open(merged_file_path, "w") as f:
        f.write(">dummy_seq\nATGCATGC\n")
    
    # When force_upload is True, the process_bins function should try to upload even if the merged file exists.
    # In this test, the simulated upload for the "fail_me" file will always fail,
    # but the process should not raise an unhandled exception.
    try:
        core.generate_merged_filtered_bins(
            min_completeness=50,
            max_contamination=10,
            min_length=1000,
            config=config
        )
    except Exception as e:
        pytest.fail(f"process_bins() raised an unexpected exception: {e}")
    
    # At this point, if the upload fails after all retries,
    # the pipeline should log an error but not crash.
    # For a more advanced test, you might capture log output to verify that the failure was noted.
    # For now, the absence of an unhandled exception is sufficient.


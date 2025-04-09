import os
import logging
import boto3
from botocore.exceptions import NoCredentialsError
import shutil
import sys

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
if not logger.handlers:
    handler = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

s3_client = boto3.client("s3")

def download_file(s3_bucket: str, s3_key: str, local_path: str) -> None:
    """
    Download a file from S3 if it does not exist locally.
    """
    if not os.path.exists(local_path):
        try:
            logger.info(f"Downloading {s3_key} from s3://{s3_bucket}/{s3_key}...")
            s3_client.download_file(s3_bucket, s3_key, local_path)
        except NoCredentialsError:
            logger.error("AWS credentials not found. Please configure your credentials.")
    else:
        logger.info(f"File {local_path} already exists. Skipping download.")

def upload_file_to_s3(local_path: str, s3_bucket: str, s3_key: str) -> None:
    """
    Upload the file at local_path to S3.
    """
    try:
        logger.info(f"Uploading {local_path} to s3://{s3_bucket}/{s3_key}")
        s3_client.upload_file(local_path, s3_bucket, s3_key)
    except NoCredentialsError:
        logger.error("AWS credentials not found. Please configure your credentials.")
    except Exception as e:
        logger.error(f"Error uploading {local_path}: {e}")

def check_tool_availability(tool: str) -> None:
    """
    Check if the given tool is available in the system PATH.
    If not, log an error and exit.
    """
    if shutil.which(tool) is None:
        error_message = (f"Error: '{tool}' is not installed or not found in your PATH. "
                         f"Please install '{tool}' and try again. For installation instructions, see: "
                         f"https://bioinf.shenwei.me/seqkit/usage/")
        import logging
        logger = logging.getLogger(__name__)
        logger.error(error_message)
        sys.exit(1)

import os
import boto3
from botocore.exceptions import NoCredentialsError

s3_client = boto3.client("s3")

def download_file(s3_bucket, s3_key, local_path):
    """
    Download a file from S3 if it does not exist locally.
    """
    if not os.path.exists(local_path):
        try:
            print(f"Downloading {s3_key} from s3://{s3_bucket}/{s3_key}...")
            s3_client.download_file(s3_bucket, s3_key, local_path)
        except NoCredentialsError:
            print("AWS credentials not found. Please configure your credentials.")
    else:
        print(f"File {local_path} already exists. Skipping download.")

def upload_file_to_s3(local_path, s3_bucket, s3_key):
    """
    Upload the file at local_path to S3.
    """
    try:
        print(f"Uploading {local_path} to s3://{s3_bucket}/{s3_key}")
        s3_client.upload_file(local_path, s3_bucket, s3_key)
    except NoCredentialsError:
        print("AWS credentials not found. Please configure your credentials.")
    except Exception as e:
        print(f"Error uploading {local_path}: {e}")


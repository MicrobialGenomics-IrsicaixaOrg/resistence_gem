import os
import boto3
import pandas as pd
import gzip
import shutil
from concurrent.futures import ThreadPoolExecutor
import subprocess
from botocore.exceptions import NoCredentialsError

# AWS S3 client
s3_client = boto3.client("s3")

############################ DOWNLOAD, FILTER, MERGE, COMPRESS AND UPLOAD CONTIGS FROM METAGENOME ASSEMBLY ###################################### 

def download_file(s3_bucket, s3_key, local_path):
    """Download a file from S3 if it does not exist locally."""
    if not os.path.exists(local_path):
        try:
            print(f"Downloading {s3_key} from s3://{s3_bucket}/{s3_key}...")
            s3_client.download_file(s3_bucket, s3_key, local_path)
        except NoCredentialsError:
            print("AWS credentials not found. Please configure your credentials.")
    else:
        print(f"File {local_path} already exists. Skipping download.")

def download_files(sample_names, assembler, s3_base_path, local_dir):
    """Download contig files from S3 in parallel."""
    os.makedirs(local_dir, exist_ok=True)
    s3_bucket = s3_base_path.split('/')[2]
    s3_prefix = "/".join(s3_base_path.split('/')[3:])
    
    with ThreadPoolExecutor() as executor:
        executor.map(lambda sample: download_file(s3_bucket, f"{s3_prefix}{assembler}-{sample}.contigs.fa.gz", 
                                                 os.path.join(local_dir, f"{assembler}-{sample}.contigs.fa.gz")), sample_names)

def filter_sequences(sample_names, assembler, local_dir, min_length=1000):
    """Filter sequences longer than min_length (default: 1000bp). [doi:10.1534/g3.118.200745]"""
    for sample in sample_names:
        input_file = os.path.join(local_dir, f"{assembler}-{sample}.contigs.fa.gz")
        output_file = os.path.join(local_dir, f"{assembler}-{sample}_filtered.fa")
        
        print(f"Filtering {input_file} (keeping sequences > {min_length} bp)...")
        subprocess.run(["seqkit", "seq", "-m", str(min_length), input_file, "-o", output_file])

def merge_filtered_files(sample_names, assembly_dir, merged_dir):
    """Merge filtered MEGAHIT and SPAdes contigs for each sample and compress the output."""
    os.makedirs(merged_dir, exist_ok=True)
    
    for sample in sample_names:
        megahit_file = os.path.join(assembly_dir, "MEGAHIT", f"MEGAHIT-{sample}_filtered.fa")
        spades_file = os.path.join(assembly_dir, "SPAdes", f"SPAdes-{sample}_filtered.fa")
        merged_file = os.path.join(merged_dir, f"{sample}_merged.fa.gz")
        
        if os.path.exists(megahit_file) and os.path.exists(spades_file):
            print(f"Merging {megahit_file} and {spades_file} into {merged_file}...")
            with gzip.open(merged_file, 'wb') as f_out:
                for file in [megahit_file, spades_file]:
                    with open(file, 'rb') as f_in:
                        shutil.copyfileobj(f_in, f_out)
        else:
            print(f"Warning: One or both files missing for {sample}, skipping merge.")

def upload_to_s3(local_path, s3_bucket, s3_key):
    """Upload the merged file to the specified S3 bucket and key."""
    try:
        print(f"Uploading {local_path} to s3://{s3_bucket}/{s3_key}")
        s3_client.upload_file(local_path, s3_bucket, s3_key)
    except NoCredentialsError:
        print("AWS credentials not found. Please configure your credentials.")

def generate_merged_filtered_contigs(samplesheet_path, min_length=1000):
    """Wrapper function to download, filter, merge, compress, and upload contigs from metagenome assembly to S3, wich will be used as input to the funcscan pipeline"""
    df = pd.read_csv(samplesheet_path)
    sample_names = df["sample"].tolist()
    
    example_path = df["short_reads_1"].iloc[0]
    s3_root_path = "/".join(example_path.split("/")[:3]) + "/"
    
    s3_megahit_path = f"{s3_root_path}Assembly/MEGAHIT/"
    s3_spades_path = f"{s3_root_path}Assembly/SPAdes/"
    s3_upload_path = f"{s3_root_path}Assembly/merged_results/"
    
    assembly_dir = "./Assembly_results/"
    merged_dir = os.path.join(assembly_dir, "merged_results")
    
    # Download and filter
    download_files(sample_names, "MEGAHIT", s3_megahit_path, os.path.join(assembly_dir, "MEGAHIT"))
    download_files(sample_names, "SPAdes", s3_spades_path, os.path.join(assembly_dir, "SPAdes"))
    
    filter_sequences(sample_names, "MEGAHIT", os.path.join(assembly_dir, "MEGAHIT"), min_length)
    filter_sequences(sample_names, "SPAdes", os.path.join(assembly_dir, "SPAdes"), min_length)
    
    merge_filtered_files(sample_names, assembly_dir, merged_dir)
    
    # Upload merged files to S3 in parallel
    s3_bucket = s3_root_path.split("/")[2]
    with ThreadPoolExecutor() as executor:
        executor.map(lambda sample: upload_to_s3(os.path.join(merged_dir, f"{sample}_merged.fa.gz"), s3_bucket, f"Assembly/merged_results/{sample}_merged.fa.gz"), sample_names)

# Run the process
generate_merged_filtered_contigs("samplesheet.csv", min_length=1000)

############################ DOWNLOAD, FILTER, MERGE, AND UPLOAD CONTIGS FROM FILTERED BINS FROM METAGENOME BINNING ######################################

def upload_to_s3(local_file_path, bucket_name, s3_file_path):
    """Upload a file to S3."""
    try:
        s3_client.upload_file(local_file_path, bucket_name, s3_file_path)
        print(f"Uploaded {local_file_path} to s3://{bucket_name}/{s3_file_path}")
    except NoCredentialsError:
        print("Credentials not available.")
    except Exception as e:
        print(f"Error uploading {local_file_path}: {e}")

# Main processing function with configurable parameters
def generate_merged_filtered_bins(min_completeness=50, max_contamination=10, min_length=1000):
    """Function to download, filter, merge, and upload contigs from bins from metagenome binning to S3, wich will be used as input to the funcscan pipeline"""
    df = pd.read_csv("samplesheet.csv")
    sample_names = df["sample"].tolist()

    # Extract S3 root path from the first short_reads_1 entry
    example_path = df["short_reads_1"].iloc[0]
    s3_root_path = "/".join(example_path.split("/")[:3]) + "/"

    print(f"Detected S3 root path: {s3_root_path}")

    # Define dynamic S3 paths
    s3_base_path_metabat = f"{s3_root_path}GenomeBinning/MetaBAT2/bins/"
    s3_base_path_maxbin = f"{s3_root_path}GenomeBinning/MaxBin2/bins/"
    s3_base_path_qc = f"{s3_root_path}GenomeBinning/QC/"
    s3_bucket = s3_root_path.split('/')[2]

    # Define output directories
    binning_dir = "./Binning_results/"
    metabat_dir = os.path.join(binning_dir, "MetaBAT")
    maxbin_dir = os.path.join(binning_dir, "MaxBin")
    merged_dir = os.path.join(binning_dir, "merged_bins")
    qc_dir = os.path.join(binning_dir, "QC_reports")

    # Ensure directories exist
    os.makedirs(metabat_dir, exist_ok=True)
    os.makedirs(maxbin_dir, exist_ok=True)
    os.makedirs(merged_dir, exist_ok=True)
    os.makedirs(qc_dir, exist_ok=True)

    # Local paths for QC reports
    busco_file_local = os.path.join(qc_dir, "busco_summary.tsv")
    quast_file_local = os.path.join(qc_dir, "quast_summary.tsv")

    # Download BUSCO & QUAST reports from S3 using boto3
    def download_file_from_s3(s3_path, local_path):
        try:
            s3_client.download_file(s3_path.split('/')[2], '/'.join(s3_path.split('/')[3:]), local_path)
            print(f"Downloaded {s3_path} to {local_path}")
        except NoCredentialsError:
            print("Credentials not available.")
        except Exception as e:
            print(f"Error downloading {s3_path}: {e}")

    print("Downloading BUSCO and QUAST reports from S3...")
    download_file_from_s3(f"{s3_base_path_qc}busco_summary.tsv", busco_file_local)
    download_file_from_s3(f"{s3_base_path_qc}quast_summary.tsv", quast_file_local)

    def filter_high_quality_bins(busco_file, quast_file, min_completeness, max_contamination):
        """Filter bins based on BUSCO completeness and QUAST contamination. By default, min_completeness=50, max_contamination=10, and min_length=1000 [doi:10.1534/g3.118.200745]"""
        print("\nFiltering high-quality bins based on BUSCO and QUAST results...")

        busco_df = pd.read_csv(busco_file, sep="\t")
        quast_df = pd.read_csv(quast_file, sep="\t")

        busco_df.rename(columns={"GenomeBin": "Bin", "%Complete (specific)": "Completeness", "%Missing (specific)": "Missing"}, inplace=True)
        quast_df.rename(columns={"Assembly": "Bin"}, inplace=True)

        if "Bin" not in busco_df.columns or "Bin" not in quast_df.columns:
            raise KeyError("Could not find a 'Bin' column in one of the files.")

        merged_df = pd.merge(busco_df, quast_df, on="Bin", how="inner")

        high_quality_bins = merged_df[
            (merged_df["Completeness"] >= min_completeness) & 
            (merged_df["Missing"] <= max_contamination)
        ]

        print(f"High-quality bins found: {len(high_quality_bins)}")
        
        # Modify the bin names to reflect the .fa.gz extension
        high_quality_bins_list = high_quality_bins["Bin"].tolist()
        high_quality_bins_gz = [f"{bin_name}.gz" for bin_name in high_quality_bins_list]
        
        return high_quality_bins_gz

    # Get high-quality bins with .fa.gz extension
    high_quality_bins = filter_high_quality_bins(busco_file_local, quast_file_local, min_completeness, max_contamination)

    # Function to filter sequences with seqkit
    def filter_sequences_with_seqkit(input_file, output_file, min_length):
        """Filter sequences based on minimum length using seqkit."""
        print(f"Filtering sequences in {input_file} with seqkit (min length: {min_length})...")
        subprocess.run(["seqkit", "seq", "-m", str(min_length), input_file, "-o", output_file])

    # Function to download, filter, save, and merge files using boto3
    def save_and_merge_files(file_list):
        merged_files = {sample: [] for sample in sample_names}  # Initialize empty lists for each sample

        for file_name in file_list:
            # Define the target subfolder based on file type
            if 'MaxBin' in file_name:
                target_subfolder = maxbin_dir
                target_s3_folder = s3_base_path_maxbin
            elif 'MetaBAT' in file_name:
                target_subfolder = metabat_dir
                target_s3_folder = s3_base_path_metabat
            else:
                print(f"Skipping unknown file type: {file_name}")
                continue

            # Ensure target folder exists (subfolder for MetaBAT or MaxBin)
            os.makedirs(target_subfolder, exist_ok=True)

            # Download file from S3 using boto3
            local_file_path = os.path.join(target_subfolder, file_name)
            if not os.path.exists(local_file_path):
                s3_url = f"{target_s3_folder}{file_name}"
                print(f"Downloading {file_name} from {s3_url}")
                try:
                    download_file_from_s3(s3_url, local_file_path)
                except Exception as e:
                    print(f"Error downloading {file_name}: {e}")
                    continue

            # Filter sequences using seqkit before merging
            filtered_file_path = local_file_path.replace(".fa.gz", "_filtered.fa.gz")
            filter_sequences_with_seqkit(local_file_path, filtered_file_path, min_length)

            # Extract the sample ID based on the file name
            for sample in sample_names:
                if sample in file_name:
                    # Add the filtered file to the list of merged files for this sample
                    merged_files[sample].append(filtered_file_path)
                    break

        # Merge files for each sample and save them in the merged_results folder
        os.makedirs(merged_dir, exist_ok=True)

        for sample_id, file_paths in merged_files.items():
            # Only merge files if there are any for this sample
            if file_paths:
                merged_file_path = os.path.join(merged_dir, f"{sample_id}_merged.fa.gz")
                
                with open(merged_file_path, 'wb') as merged_file:
                    for file_path in file_paths:
                        with open(file_path, 'rb') as f:
                            shutil.copyfileobj(f, merged_file)  # Merge files

                print(f"Merged file saved: {merged_file_path}")

                # Upload the merged file to S3
                s3_bucket = s3_root_path.split("/")[2]
                s3_file_path = f"GenomeBinning/merged_results/{sample_id}_merged.fa.gz"
                upload_to_s3(merged_file_path, s3_bucket, s3_file_path)

    # Run the function to save, merge files and upload to S3
    save_and_merge_files(high_quality_bins)

    # Parallel upload of all merged files to S3 using ThreadPoolExecutor
    with ThreadPoolExecutor() as executor:
        for sample in sample_names:
            merged_file_path = os.path.join(merged_dir, f"{sample}_merged.fa.gz")
            if os.path.exists(merged_file_path):  # Only upload if file exists
                s3_file_path = f"GenomeBinning/merged_results/{sample}_merged.fa.gz"
                executor.submit(upload_to_s3, merged_file_path, s3_bucket, s3_file_path)

generate_merged_filtered_bins(min_completeness=50, max_contamination=10, min_length=1000)

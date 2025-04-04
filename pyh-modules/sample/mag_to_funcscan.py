#!/usr/bin/env python3
"""
Process metagenomic assembly and binning results for functional analysis pipeline.

This script:
1. Downloads assembly contigs (MEGAHIT and SPAdes) from S3
2. Downloads binned contigs (MetaBAT2 and MaxBin2) from S3
3. Filters sequences based on length and quality metrics
4. Merges results from different tools
5. Uploads processed files back to S3 with 'nf_mag' subfolder

Usage:
    python mag_to_funcscan.py
"""

import os
import boto3
import pandas as pd
import gzip
import shutil
from concurrent.futures import ThreadPoolExecutor
import subprocess
from botocore.exceptions import NoCredentialsError, ClientError
import logging
from typing import List

# Configuration
CONFIG = {
    "min_contig_length": 1000,
    "min_bin_completeness": 50,
    "max_bin_contamination": 10,
    "samplesheet_path": "samplesheet.csv",
    "local_work_dir": "./processing_results",
    "nf_mag_subfolder": "nf_mag"
}

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# AWS S3 client
s3_client = boto3.client("s3")

def upload_to_s3(local_path: str, s3_bucket: str, s3_key: str) -> None:
    """Upload a file to S3 with comprehensive error handling."""
    try:
        logger.info(f"Uploading {local_path} to s3://{s3_bucket}/{s3_key}")
        s3_client.upload_file(local_path, s3_bucket, s3_key)
    except NoCredentialsError:
        logger.error("AWS credentials not found. Please configure your credentials.")
    except ClientError as e:
        logger.error(f"AWS Client Error uploading {local_path}: {e}")
    except Exception as e:
        logger.error(f"Unexpected error uploading {local_path}: {e}")

def download_file_from_s3(s3_bucket: str, s3_key: str, local_path: str) -> None:
    """Download a file from S3 with error handling."""
    try:
        if not os.path.exists(local_path):
            logger.info(f"Downloading s3://{s3_bucket}/{s3_key} to {local_path}")
            s3_client.download_file(s3_bucket, s3_key, local_path)
        else:
            logger.info(f"File {local_path} already exists. Skipping download.")
    except NoCredentialsError:
        logger.error("AWS credentials not found. Please configure your credentials.")
    except ClientError as e:
        logger.error(f"AWS Client Error downloading {s3_key}: {e}")
    except Exception as e:
        logger.error(f"Unexpected error downloading {s3_key}: {e}")

def build_s3_path(base_path: str) -> str:
    """Construct S3 path with nf_mag subfolder."""
    parts = base_path.split('/')
    bucket = parts[2]
    prefix = '/'.join(parts[3:])
    return f"s3://{bucket}/{CONFIG['nf_mag_subfolder']}/{prefix}"

def download_files(sample_names: List[str], assembler: str, s3_base_path: str, local_dir: str) -> None:
    """Download contig files from S3 in parallel."""
    os.makedirs(local_dir, exist_ok=True)
    s3_bucket = s3_base_path.split('/')[2]
    s3_prefix = "/".join(s3_base_path.split('/')[3:])
    
    with ThreadPoolExecutor() as executor:
        executor.map(
            lambda sample: download_file_from_s3(
                s3_bucket,
                f"{s3_prefix}{assembler}-{sample}.contigs.fa.gz",
                os.path.join(local_dir, f"{assembler}-{sample}.contigs.fa.gz")
            ),
            sample_names
        )

def filter_sequences(sample_names: List[str], assembler: str, local_dir: str, min_length: int) -> None:
    """Filter sequences longer than min_length."""
    for sample in sample_names:
        input_file = os.path.join(local_dir, f"{assembler}-{sample}.contigs.fa.gz")
        output_file = os.path.join(local_dir, f"{assembler}-{sample}_filtered.fa")
        
        logger.info(f"Filtering {input_file} (keeping sequences > {min_length} bp)")
        subprocess.run(
            ["seqkit", "seq", "-m", str(min_length), input_file, "-o", output_file],
            check=True
        )

def merge_filtered_files(sample_names: List[str], assembly_dir: str, merged_dir: str) -> None:
    """Merge filtered MEGAHIT and SPAdes contigs for each sample with unique sequence names."""
    os.makedirs(merged_dir, exist_ok=True)
    
    for sample in sample_names:
        megahit_file = os.path.join(assembly_dir, "MEGAHIT", f"MEGAHIT-{sample}_filtered.fa")
        spades_file = os.path.join(assembly_dir, "SPAdes", f"SPAdes-{sample}_filtered.fa")
        merged_file = os.path.join(merged_dir, f"{sample}_merged.fa.gz")
        
        if os.path.exists(megahit_file) and os.path.exists(spades_file):
            logger.info(f"Merging {megahit_file} and {spades_file} into {merged_file}")
            
            with gzip.open(merged_file, 'wt') as f_out:
                # Process MEGAHIT file
                with open(megahit_file, 'r') as f_in:
                    for line in f_in:
                        if line.startswith('>'):
                            # Add MEGAHIT prefix to sequence names
                            f_out.write(f">{sample}_MEGAHIT_{line[1:]}")
                        else:
                            f_out.write(line)
                
                # Process SPAdes file
                with open(spades_file, 'r') as f_in:
                    for line in f_in:
                        if line.startswith('>'):
                            # Add SPAdes prefix to sequence names
                            f_out.write(f">{sample}_SPAdes_{line[1:]}")
                        else:
                            f_out.write(line)
        else:
            logger.warning(f"One or both files missing for {sample}, skipping merge.")

def generate_merged_filtered_contigs(samplesheet_path: str, min_length: int) -> None:
    """Process metagenome assembly contigs."""
    logger.info("Starting assembly contigs processing")
    
    df = pd.read_csv(samplesheet_path)
    sample_names = df["sample"].tolist()
    
    example_path = df["short_reads_1"].iloc[0]
    s3_root_path = "/".join(example_path.split("/")[:3]) + "/"
    
    # Build paths with nf_mag subfolder
    s3_megahit_path = build_s3_path(f"{s3_root_path}Assembly/MEGAHIT/")
    s3_spades_path = build_s3_path(f"{s3_root_path}Assembly/SPAdes/")
    s3_upload_path = build_s3_path(f"{s3_root_path}Assembly/merged_results/")
    
    assembly_dir = os.path.join(CONFIG["local_work_dir"], "Assembly_results")
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
        executor.map(
            lambda sample: upload_to_s3(
                os.path.join(merged_dir, f"{sample}_merged.fa.gz"),
                s3_bucket,
                f"{CONFIG['nf_mag_subfolder']}/Assembly/merged_results/{sample}_merged.fa.gz"
            ),
            sample_names
        )

def filter_high_quality_bins(busco_file: str, quast_file: str, min_completeness: int, max_contamination: int) -> List[str]:
    """Filter bins based on quality metrics."""
    logger.info("Filtering high-quality bins based on BUSCO and QUAST results")
    
    busco_df = pd.read_csv(busco_file, sep="\t")
    quast_df = pd.read_csv(quast_file, sep="\t")

    busco_df.rename(columns={
        "GenomeBin": "Bin",
        "%Complete (specific)": "Completeness",
        "%Missing (specific)": "Missing"
    }, inplace=True)
    quast_df.rename(columns={"Assembly": "Bin"}, inplace=True)

    if "Bin" not in busco_df.columns or "Bin" not in quast_df.columns:
        raise KeyError("Could not find 'Bin' column in BUSCO or QUAST files")

    merged_df = pd.merge(busco_df, quast_df, on="Bin", how="inner")

    high_quality_bins = merged_df[
        (merged_df["Completeness"] >= min_completeness) & 
        (merged_df["Missing"] <= max_contamination)
    ]

    logger.info(f"Found {len(high_quality_bins)} high-quality bins")
    
    return [f"{bin_name}.gz" for bin_name in high_quality_bins["Bin"].tolist()]

def process_bin_files(file_list: List[str], sample_names: List[str], s3_base_path_metabat: str,
                     s3_base_path_maxbin: str, metabat_dir: str, maxbin_dir: str,
                     merged_dir: str, min_length: int, s3_bucket: str) -> None:
    """Process and merge bin files with unique sequence names."""
    merged_files = {sample: [] for sample in sample_names}

    for file_name in file_list:
        if 'MaxBin' in file_name:
            target_subfolder = maxbin_dir
            target_s3_folder = s3_base_path_maxbin
            bin_type = "MaxBin"
        elif 'MetaBAT' in file_name:
            target_subfolder = metabat_dir
            target_s3_folder = s3_base_path_metabat
            bin_type = "MetaBAT"
        else:
            logger.warning(f"Skipping unknown file type: {file_name}")
            continue

        os.makedirs(target_subfolder, exist_ok=True)
        local_file_path = os.path.join(target_subfolder, file_name)
        
        if not os.path.exists(local_file_path):
            s3_url = f"{target_s3_folder}{file_name}"
            logger.info(f"Downloading {file_name} from {s3_url}")
            download_file_from_s3(
                s3_url.split('/')[2],
                '/'.join(s3_url.split('/')[3:]),
                local_file_path
            )

        filtered_file_path = local_file_path.replace(".fa.gz", "_filtered.fa")
        logger.info(f"Filtering sequences in {local_file_path}")
        subprocess.run(
            ["seqkit", "seq", "-m", str(min_length), local_file_path, "-o", filtered_file_path],
            check=True
        )

        for sample in sample_names:
            if sample in file_name:
                merged_files[sample].append((filtered_file_path, bin_type))
                break

    os.makedirs(merged_dir, exist_ok=True)

    for sample_id, file_info in merged_files.items():
        if file_info:
            merged_file_path = os.path.join(merged_dir, f"{sample_id}_merged.fa.gz")
            
            with gzip.open(merged_file_path, 'wt') as merged_file:
                for file_path, bin_type in file_info:
                    with open(file_path, 'r') as f:
                        for line in f:
                            if line.startswith('>'):
                                # Add bin type prefix to sequence names
                                merged_file.write(f">{sample_id}_{bin_type}_{line[1:]}")
                            else:
                                merged_file.write(line)

            logger.info(f"Created merged file: {merged_file_path}")

            upload_to_s3(
                merged_file_path,
                s3_bucket,
                f"{CONFIG['nf_mag_subfolder']}/GenomeBinning/merged_results/{sample_id}_merged.fa.gz"
            )

def generate_merged_filtered_bins(min_completeness: int, max_contamination: int, min_length: int) -> None:
    """Process metagenome binning results."""
    logger.info("Starting bin processing")
    
    df = pd.read_csv(CONFIG["samplesheet_path"])
    sample_names = df["sample"].tolist()
    example_path = df["short_reads_1"].iloc[0]
    s3_root_path = "/".join(example_path.split("/")[:3]) + "/"

    # Build paths with nf_mag subfolder
    s3_base_path_metabat = build_s3_path(f"{s3_root_path}GenomeBinning/MetaBAT2/bins/")
    s3_base_path_maxbin = build_s3_path(f"{s3_root_path}GenomeBinning/MaxBin2/bins/")
    s3_base_path_qc = build_s3_path(f"{s3_root_path}GenomeBinning/QC/")
    s3_bucket = s3_root_path.split('/')[2]

    binning_dir = os.path.join(CONFIG["local_work_dir"], "Binning_results")
    metabat_dir = os.path.join(binning_dir, "MetaBAT")
    maxbin_dir = os.path.join(binning_dir, "MaxBin")
    merged_dir = os.path.join(binning_dir, "merged_bins")
    qc_dir = os.path.join(binning_dir, "QC_reports")

    os.makedirs(metabat_dir, exist_ok=True)
    os.makedirs(maxbin_dir, exist_ok=True)
    os.makedirs(merged_dir, exist_ok=True)
    os.makedirs(qc_dir, exist_ok=True)

    busco_file_local = os.path.join(qc_dir, "busco_summary.tsv")
    quast_file_local = os.path.join(qc_dir, "quast_summary.tsv")

    logger.info("Downloading QC reports")
    download_file_from_s3(
        s3_base_path_qc.split('/')[2],
        f"{CONFIG['nf_mag_subfolder']}/GenomeBinning/QC/busco_summary.tsv",
        busco_file_local
    )
    download_file_from_s3(
        s3_base_path_qc.split('/')[2],
        f"{CONFIG['nf_mag_subfolder']}/GenomeBinning/QC/quast_summary.tsv",
        quast_file_local
    )

    high_quality_bins = filter_high_quality_bins(
        busco_file_local,
        quast_file_local,
        min_completeness,
        max_contamination
    )

    process_bin_files(
        high_quality_bins,
        sample_names,
        s3_base_path_metabat,
        s3_base_path_maxbin,
        metabat_dir,
        maxbin_dir,
        merged_dir,
        min_length,
        s3_bucket
    )

def cleanup(local_dirs_to_remove: List[str]) -> None:
    """Remove local directories after processing."""
    for dir_path in local_dirs_to_remove:
        if os.path.exists(dir_path):
            shutil.rmtree(dir_path)
            logger.info(f"Removed local directory: {dir_path}")

def main() -> None:
    """Main execution function."""
    try:
        # Process assembly contigs
        generate_merged_filtered_contigs(
            CONFIG["samplesheet_path"],
            CONFIG["min_contig_length"]
        )
        
        # Process binned contigs
        generate_merged_filtered_bins(
            min_completeness=CONFIG["min_bin_completeness"],
            max_contamination=CONFIG["max_bin_contamination"],
            min_length=CONFIG["min_contig_length"]
        )
        
        # Cleanup (optional - comment out if you want to keep intermediate files)
        # cleanup([
        #     os.path.join(CONFIG["local_work_dir"], "Assembly_results"),
        #     os.path.join(CONFIG["local_work_dir"], "Binning_results")
        # ])
        
    except Exception as e:
        logger.error(f"Script failed: {e}", exc_info=True)
        raise

if __name__ == "__main__":
    main()

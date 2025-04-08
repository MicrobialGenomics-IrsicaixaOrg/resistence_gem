import os
import subprocess
import gzip
import shutil
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
import boto3
from botocore.exceptions import NoCredentialsError
from .helpers import download_file, upload_file_to_s3
from .constants import DEFAULT_CONFIG

# Initialize S3 client
s3_client = boto3.client("s3")


def build_s3_path(base_path: str) -> str:
    """Construct S3 path with nf_mag subfolder."""
    parts = base_path.split('/')
    bucket = parts[2]
    prefix = '/'.join(parts[3:])
    return f"s3://{bucket}/{DEFAULT_CONFIG['nf_mag_subfolder']}/{prefix}"


def download_files(sample_names, assembler, s3_base_path, local_dir):
    """
    Download contig files from S3 in parallel.
    """
    os.makedirs(local_dir, exist_ok=True)
    s3_bucket = s3_base_path.split('/')[2]
    s3_prefix = "/".join(s3_base_path.split('/')[3:])

    def download_sample(sample):
        key = f"{s3_prefix}{assembler}-{sample}.contigs.fa.gz"
        local_path = os.path.join(local_dir, f"{assembler}-{sample}.contigs.fa.gz")
        download_file(s3_bucket, key, local_path)

    with ThreadPoolExecutor() as executor:
        executor.map(download_sample, sample_names)


def filter_sequences(sample_names, assembler, local_dir, min_length=1000):
    """
    Filter sequences longer than min_length using seqkit.
    """
    for sample in sample_names:
        input_file = os.path.join(local_dir, f"{assembler}-{sample}.contigs.fa.gz")
        output_file = os.path.join(local_dir, f"{assembler}-{sample}_filtered.fa")
        print(f"Filtering {input_file} (keeping sequences > {min_length} bp)...")
        subprocess.run(["seqkit", "seq", "-m", str(min_length), input_file, "-o", output_file], check=True)


def merge_filtered_files(sample_names, assembly_dir, merged_dir):
    """
    Merge filtered contig files for each sample and compress the output.
    This version prefixes sequence headers with the assembler name so that they become unique.
    """
    os.makedirs(merged_dir, exist_ok=True)

    for sample in sample_names:
        megahit_file = os.path.join(assembly_dir, "MEGAHIT", f"MEGAHIT-{sample}_filtered.fa")
        spades_file = os.path.join(assembly_dir, "SPAdes", f"SPAdes-{sample}_filtered.fa")
        merged_file = os.path.join(merged_dir, f"{sample}_merged.fa.gz")

        if os.path.exists(megahit_file) and os.path.exists(spades_file):
            print(f"Merging {megahit_file} and {spades_file} into {merged_file}...")
            with gzip.open(merged_file, 'wt') as f_out:
                # Process MEGAHIT file: add "MEGAHIT_" prefix to sequence headers
                with open(megahit_file, 'r') as f_in:
                    for line in f_in:
                        if line.startswith('>'):
                            new_header = f">{sample}_MEGAHIT_{line[1:].strip()}\n"
                            f_out.write(new_header)
                        else:
                            f_out.write(line)
                # Process SPAdes file: add "SPAdes_" prefix to sequence headers
                with open(spades_file, 'r') as f_in:
                    for line in f_in:
                        if line.startswith('>'):
                            new_header = f">{sample}_SPAdes_{line[1:].strip()}\n"
                            f_out.write(new_header)
                        else:
                            f_out.write(line)
        else:
            print(f"Warning: One or both files missing for {sample}, skipping merge.")


def generate_merged_filtered_contigs(samplesheet_path, min_length=1000, config=DEFAULT_CONFIG):
    """
    Process contigs: download, filter, merge and upload.
    """
    samplesheet_path = config.get("samplesheet_path", "samplesheet.csv")
    df = pd.read_csv(samplesheet_path)
    sample_names = df["sample"].tolist()
    example_path = df["short_reads_1"].iloc[0]
    s3_root_path = "/".join(example_path.split("/")[:3]) + "/"

    s3_megahit_path = build_s3_path(f"{s3_root_path}Assembly/MEGAHIT/")
    s3_spades_path = build_s3_path(f"{s3_root_path}Assembly/SPAdes/")

    local_work_dir = config.get("local_work_dir", os.getcwd())
    assembly_dir = os.path.join(local_work_dir, "Assembly_results")
    merged_dir_out = os.path.join(assembly_dir, "merged_results")

    # Download files for both assemblers
    download_files(sample_names, "MEGAHIT", s3_megahit_path, os.path.join(assembly_dir, "MEGAHIT"))
    download_files(sample_names, "SPAdes", s3_spades_path, os.path.join(assembly_dir, "SPAdes"))

    # Filter sequences
    filter_sequences(sample_names, "MEGAHIT", os.path.join(assembly_dir, "MEGAHIT"), min_length)
    filter_sequences(sample_names, "SPAdes", os.path.join(assembly_dir, "SPAdes"), min_length)

    # Merge filtered files
    merge_filtered_files(sample_names, assembly_dir, merged_dir_out)

    # Upload merged files to S3
    s3_bucket = s3_root_path.split("/")[2]

    def upload(sample):
        local_file = os.path.join(merged_dir_out, f"{sample}_merged.fa.gz")
        s3_key = f"{DEFAULT_CONFIG['nf_mag_subfolder']}/Assembly/merged_results/{sample}_merged.fa.gz"
        upload_file_to_s3(local_file, s3_bucket, s3_key)

    with ThreadPoolExecutor() as executor:
        executor.map(upload, sample_names)


def generate_merged_filtered_bins(min_completeness=50, max_contamination=10, min_length=1000, config=DEFAULT_CONFIG):
    """
    Process bins: download QC reports, filter bins, download and merge bin files, then upload merged bins.
    """
    samplesheet_path = config.get("samplesheet_path", "samplesheet.csv")
    df = pd.read_csv(samplesheet_path)
    sample_names = df["sample"].tolist()
    example_path = df["short_reads_1"].iloc[0]
    s3_root_path = "/".join(example_path.split("/")[:3]) + "/"

    s3_base_path_metabat = build_s3_path(f"{s3_root_path}GenomeBinning/MetaBAT2/bins/")
    s3_base_path_maxbin = build_s3_path(f"{s3_root_path}GenomeBinning/MaxBin2/bins/")
    s3_base_path_qc = build_s3_path(f"{s3_root_path}GenomeBinning/QC/")
    s3_bucket = s3_root_path.split('/')[2]

    # Define local directories
    local_work_dir = config.get("local_work_dir", os.getcwd())
    binning_dir = os.path.join(local_work_dir, "Binning_results")
    metabat_dir = os.path.join(binning_dir, "MetaBAT")
    maxbin_dir = os.path.join(binning_dir, "MaxBin")
    merged_dir_bins = os.path.join(binning_dir, "merged_bins")
    qc_dir = os.path.join(binning_dir, "QC_reports")

    for d in [metabat_dir, maxbin_dir, merged_dir_bins, qc_dir]:
        os.makedirs(d, exist_ok=True)

    busco_file_local = os.path.join(qc_dir, "busco_summary.tsv")
    quast_file_local = os.path.join(qc_dir, "quast_summary.tsv")

    # Download QC reports
    def download_qc(s3_path, local_path):
        parts = s3_path.split('/')
        bucket = parts[2]
        key = "/".join(parts[3:])
        try:
            s3_client.download_file(bucket, key, local_path)
            print(f"Downloaded {s3_path} to {local_path}")
        except NoCredentialsError:
            print("AWS credentials not found.")
        except Exception as e:
            print(f"Error downloading {s3_path}: {e}")

    download_qc(f"{s3_base_path_qc}busco_summary.tsv", busco_file_local)
    download_qc(f"{s3_base_path_qc}quast_summary.tsv", quast_file_local)

    # Filter high quality bins based on BUSCO and QUAST data
    print("Filtering high-quality bins...")
    busco_df = pd.read_csv(busco_file_local, sep="\t")
    quast_df = pd.read_csv(quast_file_local, sep="\t")

    busco_df.rename(columns={"GenomeBin": "Bin", "%Complete (specific)": "Completeness", "%Missing (specific)": "Missing"}, inplace=True)
    quast_df.rename(columns={"Assembly": "Bin"}, inplace=True)

    merged_df = pd.merge(busco_df, quast_df, on="Bin", how="inner")

    high_quality_bins = merged_df[
        (merged_df["Completeness"] >= min_completeness) &
        (merged_df["Missing"] <= max_contamination)
    ]
    high_quality_bins_list = high_quality_bins["Bin"].tolist()
    high_quality_bins_gz = [f"{bin_name}.gz" for bin_name in high_quality_bins_list]

    # Process bin files: download, filter and merge
    def process_bin_files(file_list, sample_names, s3_base_path_metabat, s3_base_path_maxbin,
                          metabat_dir, maxbin_dir, merged_dir_bins, min_length, s3_bucket):
        """
        Process bin files: download, filter, and merge files per sample.
        Each filtered file is annotated with its bin type to create unique headers in the final merged file.
        """
        import subprocess, os, shutil, gzip
        from concurrent.futures import ThreadPoolExecutor

        # Create a dictionary to accumulate filtered file paths along with their bin type per sample.
        merged_files = {sample: [] for sample in sample_names}

        def process_bin(file_name):
            if "MaxBin" in file_name:
                target_dir = maxbin_dir
                target_s3_folder = s3_base_path_maxbin
                bin_type = "MaxBin"
            elif "MetaBAT" in file_name:
                target_dir = metabat_dir
                target_s3_folder = s3_base_path_metabat
                bin_type = "MetaBAT"
            else:
                print(f"Skipping unknown file type: {file_name}")
                return None, None

            os.makedirs(target_dir, exist_ok=True)
            local_file = os.path.join(target_dir, file_name)
            if not os.path.exists(local_file):
                s3_url = f"{target_s3_folder}{file_name}"
                parts = s3_url.split('/')
                bucket = parts[2]
                key = "/".join(parts[3:])
                try:
                    s3_client.download_file(bucket, key, local_file)
                    print(f"Downloaded {file_name} from {s3_url}")
                except Exception as e:
                    print(f"Error downloading {file_name}: {e}")
                    return None, None

            # Create a filtered version of the file
            filtered_file = local_file.replace(".fa.gz", "_filtered.fa.gz")
            print(f"Filtering bin file {local_file} (min length: {min_length})...")
            subprocess.run(["seqkit", "seq", "-m", str(min_length), local_file, "-o", filtered_file],
                           check=True)
            return filtered_file, bin_type

        # Process each bin file in the provided list.
        for file_name in file_list:
            filtered, bin_type = process_bin(file_name)
            if filtered:
                for sample in sample_names:
                    if sample in file_name:
                        merged_files[sample].append((filtered, bin_type))
                        break

        os.makedirs(merged_dir_bins, exist_ok=True)

        # For each sample, merge its filtered files into one merged file with modified headers.
        for sample, file_info in merged_files.items():
            if file_info:
                merged_file_path = os.path.join(merged_dir_bins, f"{sample}_merged.fa.gz")
                # Open merged file in text mode to process headers, output as gzipped file.
                with gzip.open(merged_file_path, 'wt') as wfd:
                    for file_path, bin_type in file_info:
                        # Open the filtered file using gzip in text mode.
                        with gzip.open(file_path, 'rt') as f_in:
                            for line in f_in:
                                if line.startswith('>'):
                                    new_header = f">{sample}_{bin_type}_{line[1:].strip()}\n"
                                    wfd.write(new_header)
                                else:
                                    wfd.write(line)
                print(f"Merged bin file saved: {merged_file_path}")

        # Define a helper function for parallel uploading.
        def upload_sample(sample):
            merged_file_path = os.path.join(merged_dir_bins, f"{sample}_merged.fa.gz")
            if os.path.exists(merged_file_path):
                s3_key = f"{DEFAULT_CONFIG['nf_mag_subfolder']}/GenomeBinning/merged_results/{sample}_merged.fa.gz"
                upload_file_to_s3(merged_file_path, s3_bucket, s3_key)

        # Use ThreadPoolExecutor to parallelize uploading for each sample.
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor() as executor:
            executor.map(upload_sample, sample_names)

    process_bin_files(
        high_quality_bins_gz,
        sample_names,
        s3_base_path_metabat,
        s3_base_path_maxbin,
        metabat_dir,
        maxbin_dir,
        merged_dir_bins,
        min_length,
        s3_bucket
    )


def main():
    """
    Command-line interface for the library.
    """
    import argparse
    parser = argparse.ArgumentParser(
        description="Process metagenome assemblies: contigs and bins merging pipeline"
    )
    subparsers = parser.add_subparsers(dest="command", help="Sub-command help")

    contigs_parser = subparsers.add_parser("contigs", help="Process contigs")
    contigs_parser.add_argument("samplesheet", help="Path to the samplesheet CSV file")
    contigs_parser.add_argument("--min_length", type=int, default=1000, help="Minimum sequence length for filtering")

    bins_parser = subparsers.add_parser("bins", help="Process bins")
    bins_parser.add_argument("--min_completeness", type=int, default=50, help="Minimum BUSCO completeness")
    bins_parser.add_argument("--max_contamination", type=int, default=10, help="Maximum QUAST contamination")
    bins_parser.add_argument("--min_length", type=int, default=1000, help="Minimum sequence length for filtering")

    args = parser.parse_args()

    if args.command == "contigs":
        generate_merged_filtered_contigs(args.samplesheet, min_length=args.min_length)
    elif args.command == "bins":
        generate_merged_filtered_bins(
            min_completeness=args.min_completeness,
            max_contamination=args.max_contamination,
            min_length=args.min_length
        )
    else:
        parser.print_help()


if __name__ == "__main__":
    main()


class MAGProcessor:
    def __init__(self, config=None):
        # Merge the provided config with the default configuration.
        self.config = DEFAULT_CONFIG.copy()
        if config:
            self.config.update(config)

    def process_assemblies(self):
        """
        Process assembly contigs: download, filter, merge, and upload.
        """
        from .core import generate_merged_filtered_contigs
        generate_merged_filtered_contigs(
            self.config["samplesheet_path"],
            min_length=self.config.get("min_contig_length", 1000),
            config=self.config
        )

    def process_bins(self):
        """
        Process bins: download, filter, merge, and upload.
        """
        from .core import generate_merged_filtered_bins
        generate_merged_filtered_bins(
            min_completeness=self.config.get("min_bin_completeness", 50),
            max_contamination=self.config.get("max_bin_contamination", 10),
            min_length=self.config.get("min_contig_length", 1000),
            config=self.config
        )


pyh-modules
===========

A Python library for processing metagenomic assembly and binning results from S3.

Features
--------
- Download assembly contigs (MEGAHIT/SPAdes) and bins (MetaBAT2/MaxBin2) from S3.
- Filter sequences by length and quality metrics.
- Merge results from multiple tools.
- Upload processed files back to S3.

Installation
------------
.. code-block:: bash

    # Clone the repository
    git clone https://github.com/MicrobialGenomics-IrsicaixaOrg/resistence_gem.git
    cd resistence_gem/pyh-modules

    # Create and activate a virtual environment (recommended)
    python -m venv venv
    source venv/bin/activate   # Linux/Mac
    # OR
    venv\Scripts\activate      # Windows

    # Install the package in editable (development) mode
    pip install -e .

Configuration
-------------
### AWS Credentials Setup

1. **Using AWS CLI (Recommended):**
   - Install the AWS CLI:  
     https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html
   - Configure your credentials by running:

     .. code-block:: bash

        aws configure

   - Provide your AWS Access Key ID, AWS Secret Access Key, default region (e.g., `eu-west-1`), and output format (`json`).

   - **Verification:**
     
     .. code-block:: bash

        aws sts get-caller-identity

2. **Alternative Methods:**
   - **Environment Variables:**

     .. code-block:: bash

        export AWS_ACCESS_KEY_ID="your_access_key"
        export AWS_SECRET_ACCESS_KEY="your_secret_key"
        export AWS_DEFAULT_REGION="your_region"

   - **Shared Credentials File:**
     Edit `~/.aws/credentials` (Linux/Mac) or `%UserProfile%\.aws\credentials` (Windows):

     .. code-block:: ini

        [default]
        aws_access_key_id = YOUR_ACCESS_KEY
        aws_secret_access_key = YOUR_SECRET_KEY
        region = your_region

Troubleshooting:
- Ensure your credentials grant S3 read/write access (e.g., `AmazonS3FullAccess`).
- For CLI errors, refer to:
  https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-troubleshooting.html

### Samplesheet Format

Prepare a CSV file with the following columns:

.. code-block:: text

    sample,group,short_reads_1,short_reads_2
    sample1,0,s3://bucket/path/sample1_R1.fastq.gz,s3://bucket/path/sample1_R2.fastq.gz
    sample2,1,s3://bucket/path/sample2_R1.fastq.gz,s3://bucket/path/sample2_R2.fastq.gz

You can place the samplesheet anywhere (for example, in the repository root or in a dedicated data folder) and provide its path when running the pipeline.

### Default Configurations

Default settings are defined in ``sample/constants.py``. Key parameters include:

.. code-block:: python

    DEFAULT_CONFIG = {
        "min_contig_length": 1000,      # Minimum contig length (bp)
        "min_bin_completeness": 50,     # Minimum bin completeness (%)
        "max_bin_contamination": 10,    # Maximum bin contamination (%)
        "samplesheet_path": "samplesheet.csv",  # Path to the samplesheet
        "local_work_dir": "./processing_results",  # Directory for processing results
        "nf_mag_subfolder": "nf_mag"    # S3 subfolder for processed files
    }

Usage
-----
The pipeline is accessible via a command-line interface:

.. code-block:: bash

    # Process contigs (assemblies)
    resistence_gem contigs samplesheet.csv --min_length 1000

    # Process bins
    resistence_gem bins --min_completeness 50 --max_contamination 10 --min_length 1000

Alternatively, you can use the library programmatically:

.. code-block:: python

    from sample import MAGProcessor

    # Using default configuration
    processor = MAGProcessor()

    # Or, using a custom configuration
    custom_config = {
        "min_contig_length": 1500,
        "samplesheet_path": "data/my_samples.csv"
    }
    processor = MAGProcessor(config=custom_config)

    # Run the pipeline
    processor.process_assemblies()  # Process contigs
    processor.process_bins()        # Process bins

File Structure
--------------
.. code-block:: text

    resistence_gem/
    ├── pyh-modules/
    │   ├── docs/                     # Documentation (Sphinx)
    │   ├── sample/                   # Library code
    │   │   ├── __init__.py           # Package initialization
    │   │   ├── constants.py          # Default configuration values
    │   │   ├── core.py               # Main processing logic
    │   │   ├── helpers.py            # AWS S3 utilities
    │   │   └── exceptions.py         # Custom exceptions
    │   ├── tests/                    # Unit and integration tests
    │   ├── samplesheet.csv           # Example samplesheet (optional)
    │   ├── requirements.txt          # Package dependencies
    │   └── setup.py                  # Package setup script

Testing
-------
Run tests using pytest. From the repository root, execute:

.. code-block:: bash

    pytest -v

Support
-------
### Reporting Issues
1. **GitHub Issues (Preferred):**
   - Navigate to: https://github.com/MicrobialGenomics-IrsicaixaOrg/resistence_gem/issues
   - Click "New Issue" and select the appropriate template.
   - Include error logs, steps to reproduce, and the version of ``pyh-modules`` 
     (check with: `python -c "from sample import __version__; print(__version__)"`).

2. **Direct Contact:**
   - Email: ocareta@irsicaixa.es
   - Subject: `[resistence_gem/pyh-modules] Issue Description`

### Access Notes
- This repository is private. Contact the maintainers for access to:
  - Code: https://github.com/MicrobialGenomics-IrsicaixaOrg/resistence_gem/tree/main/pyh-modules
  - Internal documentation (if hosted separately)


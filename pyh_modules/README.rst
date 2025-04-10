pyh_modules
===========

A Python library for processing metagenomic assembly and binning results from S3.

Features
--------
- Download assembly contigs (MEGAHIT/SPAdes) and bins (MetaBAT2/MaxBin2) from S3.
- Filter sequences by length and quality metrics.
- Merge results from multiple tools.
- Support for resumable uploads: if an upload fails (e.g., due to temporary network or SSL issues),
  the pipeline will log the error and retry the upload.
- Configurable re-upload: use the *force_upload* parameter to force re-upload even when local
  merged files exist.

Installation
------------
.. code-block:: bash

    # Clone the repository
    git clone https://github.com/MicrobialGenomics-IrsicaixaOrg/resistence_gem.git
    cd resistence_gem/pyh_modules

    # Create and activate a virtual environment (recommended)
    python -m venv venv
    source venv/bin/activate   # Linux/Mac
    # OR
    venv\Scripts\activate      # Windows

    # Install the package in editable (development) mode
    pip install -e .
    

**Note:** Ensure that **seqkit** is installed and available on your system's PATH. Check out the installation instructions at:
https://bioinf.shenwei.me/seqkit/usage/

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

The default configuration in ``pyh_modules/constants.py`` includes parameters such as:

.. code-block:: python

    DEFAULT_CONFIG = {
        "min_contig_length": 1000,
        "min_bin_completeness": 50,
        "max_bin_contamination": 10,
        "samplesheet_path": "samplesheet.csv",
        "local_work_dir": "./processing_results",
        "nf_mag_subfolder": "nf_mag",
        "force_upload": False  # If True, forces re-upload of merged files, even if they exist locally.
    }

By enabling *force_upload*, both the assembly and bin processing pipelines will attempt to re-upload files (using retries) even if the merged files are already present locally. This allows you to overwrite files on S3 in case of previous upload failures or if you intentionally want to refresh the remote copy.

To force re-upload (for example, if you want to overwrite existing files on S3), set the configuration:

.. code-block:: python

    custom_config = {
        "samplesheet_path": "/path/to/your/samplesheet.csv",
        "local_work_dir": "/path/to/your/results",
        "force_upload": True
    }
    processor = MAGProcessor(config=custom_config)
    processor.process_bins()

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

    python3
    
    from pyh_modules import MAGProcessor

    # Using default configuration
    processor = MAGProcessor()

    # Or, using a custom configuration
        custom_config = {
        "samplesheet_path": "/home/user/Documents/data/samplesheet.csv",
        "local_work_dir": "/home/user/Documents/data/results",
        "force_upload": True  # Force re-upload if needed
    }
 
    processor = MAGProcessor(config=custom_config)

    # Run the pipeline
    processor.process_assemblies()  # Process contigs
    processor.process_bins()        # Process bins

File Structure
--------------
.. code-block:: text

    resistence_gem/
    ├── pyh_modules/                  # Root directory of the package
    ├── docs/                       # Documentation files (Sphinx configuration and index)
    │   ├── conf.py                 # Sphinx configuration file
    │   └── index.rst               # Main documentation index
    ├── README.rst                  # Project README with installation, configuration, and usage instructions
    ├── requirements.txt            # List of Python package dependencies
    ├── pyh_modules/                     # Main package code
    │   ├── constants.py            # Default configuration settings
    │   ├── core.py                 # Core processing logic (e.g., downloads, filtering, merging, uploads)
    │   ├── helpers.py              # Helper functions (e.g., S3 operations)
    │   ├── __init__.py             # Package initialization file (exports package components)
    ├── setup.py                    # Setup script for packaging and installation
    └── tests/                      # Test suite for unit and integration tests
        ├── conftest.py             # Pytest configuration and fixtures
        ├── test_advanced.py        # Advanced integration tests
        └── test_basic.py           # Basic unit tests

        
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


from setuptools import setup, find_packages

# Use robust file reading for the README
try:
    with open("README.rst", encoding="utf-8") as fh:
        long_description = fh.read()
except Exception as e:
    long_description = ""
    print(f"Warning: Could not read README.rst: {e}")

setup(
    name="pyh_modules",
    version="0.1.0",
    # Updated to search in the renamed package folder (update directory name accordingly)
    packages=find_packages(where="pyh_modules"),
    package_dir={"": "pyh_modules"},
    install_requires=[
        "boto3>=1.20.0",
        "pandas>=1.3.0"
    ],
    entry_points={
        "console_scripts": [
            "resistence_gem=pyh_modules.core:main",
        ],
    },
    author="OriolGEM",
    author_email="ocareta@irsicaixa.es",
    description="A Python library for processing metagenome assemblies: downloading, filtering, merging, and uploading contigs/bins.",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/MicrobialGenomics-IrsicaixaOrg/resistence_gem",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)


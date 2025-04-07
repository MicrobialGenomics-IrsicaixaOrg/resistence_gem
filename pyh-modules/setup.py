from setuptools import setup, find_packages

setup(
    name="pyh_modules",
    version="0.1.0",
    packages=find_packages(where="sample"),
    package_dir={"": "sample"},
    install_requires=[
        "boto3",
        "pandas"
    ],
    entry_points={
        "console_scripts": [
            "resistence_gem=core:main",
        ],
    },
    author="OriolGEM",
    author_email="ocareta@irsicaixa.es",
    description="A Python library for processing metagenome assemblies: downloading, filtering, merging, and uploading contigs/bins.",
    long_description=open("README.rst").read(),
    long_description_content_type="text/x-rst",
    url="https://github.com/MicrobialGenomics-IrsicaixaOrg/resistence_gem",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)


from setuptools import setup, find_packages

setup(
    name="HyphaeS",
    version="1.0.0",
    author="Shrihari Kamalan Kumaraguruparan",
    author_email="kkshrihari@gmail.com",
    description="Shotgun soil metagenomics pipeline for low-abundance fungi",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourname/HyphaeS",
    license="MIT",
    packages=find_packages(),
    package_data={
        "hyphaes": [
            "templates/*.yaml",
            "templates/*.tsv",
            "workflow/Snakefile",
            "workflow/rules/*.smk",
            "workflow/envs/*.yaml",
            "workflow/schemas/*.yaml",
            "workflow/scripts/*.py",
        ]
    },
    install_requires=[
        "click>=8.0",
        "snakemake>=7.0.0",
        "pandas>=1.3",
        "numpy>=1.22.4",
        "matplotlib>=3.4",
        "pyyaml>=5.4",
        "jsonschema>=4.0",
        "pulp==2.7.0",
    ],
    entry_points={
        "console_scripts": [
            "HyphaeS=hyphaes.cli:main",
        ]
    },
    python_requires=">=3.8",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)

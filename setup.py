import sys
from setuptools import setup
from warnings import warn

if sys.version_info.major != 3:
    raise RuntimeError("Palantir requires Python 3")
if sys.version_info.minor < 6:
    warn("Analysis methods were developed using Python 3.6")

# get version
with open("src/palantir/version.py") as f:
    exec(f.read())
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="palantir",
    version=__version__,  # read in from the exec of version.py; ignore error
    description="Palantir for modeling continuous cell state and cell fate choices in single cell data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dpeerlab/palantir",
    author=__author__,
    author_email=__author_email__,
    package_dir={"": "src"},
    packages=["palantir"],
    install_requires=[
        "numpy==1.24.2",
        "pandas==1.5.3",
        "scipy==1.10.1",
        "networkx>=2.1",
        "scikit-learn>=0.24",
        "joblib",
        "fcsparser>=0.2.7",
        "leidenalg>=0.9.1",
        "Cython",
        "cmake",
        "matplotlib==3.7.3",
        "tzlocal",
        "anndata>=0.8.0",
        "scanpy>=1.6.0",
        "mellon>=1.3.0",
        "pygam",
        "pybind11"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Operating System :: POSIX :: Linux",
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    python_requires=">=3.6",
)

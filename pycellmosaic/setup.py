from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as f:
    requirements = f.read().splitlines()

setup(
    name="pycellmosaic",
    version="0.1.0",
    author="Anas Saedwi",
    author_email="author@example.com",
    description="Mosaic views for multi-omics single-cell data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/anas-saedwi/pycellmosaic",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.9",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "pycellmosaic=pycellmosaic.cli:app",
        ],
    },
)

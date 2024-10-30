#!/usr/bin/env python3

import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

d = {}
with open("temsomap/_version_.py") as f:
    exec(f.read(), d)

setuptools.setup(
    name="TemSOMap",
    version=d["__version__"],
    author="Xinhai Pan",
    include_package_data=True,
    author_email="xpan78@gatech.edu",
    description="Temporal-Spatial-Omics Mapping of single cells.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_namespace_packages(),
    classifiers=["Programming Language :: Python :: 3.11", "Operating System :: Unix",],
    python_requires=">=3.8.5",
    install_requires=[
        "pip",
        "torch",
        "pandas",
        "numpy",
        "scipy",
        "matplotlib",
        "seaborn",
        "scanpy",
        "scikit-learn",
        "tqdm",
    ],
)
#!/usr/bin/env python3
import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="worm",
    version="0.0.2",
    author="Fabian Dechent",
    author_email="dechent.fabian@gmail.com",
    package_dir={"": "worm"},
    packages=setuptools.find_packages(where="worm"),
    python_requires=">=3.8",
)

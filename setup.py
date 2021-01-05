# anonymize_bam.py - de-identify sequencing data

import setuptools

# Set __version__ for the project.
exec(open("./anonymizeBAM/version.py").read())

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="anonymizeBAM",
    version=__version__,
    author="Christoph Ziegenhain",
    author_email="christoph.ziegenhain@ki.se",
    description="de-identification sequencing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sandberg-lab/dataprivacy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3.0"
    ],
    python_requires='>=3.6',
    license='GPL-3.0-or-later'
)

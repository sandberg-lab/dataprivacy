
## anonymizeBAM.py: de-identification of sequencing reads

**anonymizeBAM.py** is a tool that can sanitize sequencing reads stored in BAM file format to protect the privacy and genetic information of donor individuals.

## Installation

anonymizeBAM.py is available through PyPI. To install, type the following command line, and add `-U` for upgrading:
`pip install -U vireoSNP`

Alternatively, you can install from this GitHub repository for the latest version:
`pip install -U git+https://github.com/sandberg-lab/dataprivacy`

Add `--user` if you don't have write permissions in your default python folder.

## Usage

anonymizeBAM.py requires only an aligned .bam file and the reference genome in fasta format.
Your fasta file should be indexed (`samtools faidx`).
The .bam file should be coordinate sorted and indexed, however `anonymizeBAM.py` will try to do this for you if not.


    usage: anonymizeBAM.py [-h] [--bam FILENAME] [--out FILENAME] [--fa FILENAME]
                            [--p P] [--strict] [--keepunmapped]
    
    optional arguments:
      -h, --help      show this help message and exit
      --bam FILENAME  Path to input BAM file
      --out FILENAME  Path to output bam file
      --fa FILENAME   Path to genome reference fasta
      --p P           Number of processes to use
      --strict        Strict: also sanitize mapping score & auxiliary tags (eg. AS / NH).
      --keepunmapped  Keep ummapped reads in output bam file.

## Description

xyz

## Reference

[link to follow]

#!/usr/bin/env python3
"""Takes a VCF file and converts it to hdf5 format using scikit-allel.

Usage: ./vcf_to_h5.py [vcf_file]

Where:
[vcf_file] is the full filepath to the VCF we want to convert
"""

import sys
import os
import numpy as np
import scipy
import pandas
import h5py
import allel

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)


def main(vcf_file):
    """Driver function."""
    tmp_filename = os.path.basename(vcf_file)
    h5_file = ''.join([os.path.dirname(vcf_file), '/', os.path.splitext(tmp_filename)[0], '.h5'])
    # Convert VCF to hdf5
    allel.vcf_to_hdf5(vcf_file, h5_file, fields='*', overwrite=True)


main(sys.argv[1]) # Run program

#!/usr/bin/env python
"""
    fnu_to_flam.py

    Converts an fnu spectrum (wavelength in Angstroms and flux in ergs/sec/cm**2/Hz) 
    into an flam spectrum (wavelength in Angstroms and flux in ergs/sec/cm**2/Angstrom).

    The spectrum file is assumed to be in CSV format, with columns "wave" and "flux".

    Examples:

    fnu_to_flam.py --help

    fnu_to_flam.py --inputFileName SSSJ0005-0127_sum1.txt.csv --outputFileName SSSJ0005-0127_sum1_flam.txt.csv --verbose 1

    """


# Initial setup...
import numpy as np
import pysynphot as S
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn
import sys
import os

#--------------------------------------------------------------------------

# Main code.

def main():

    import argparse
    import warnings
    from astropy.utils.exceptions import AstropyWarning

    # Ignore Astropy warnings...
    warnings.simplefilter('ignore', category=AstropyWarning)
    # Ignore FutureWarnings...
    warnings.simplefilter('ignore', category=FutureWarning)
    
    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputFileName', help='name of the input CSV file')
    parser.add_argument('--outputFileName', help='name of the output CSV file')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    status = fnu_to_flam(args)

    return status

#--------------------------------------------------------------------------

def fnu_to_flam(args):

    #  Extract the input spectrum filename...
    inputFileName = args.inputFileName
    if os.path.isfile(inputFileName)==False:
        print """inputFileName %s does not exist...""" % (inputFileName)
        print 'Returning with error code 1 now...'
        return 1
    if args.verbose > 0:
        print 'inputFileName: ', inputFileName

    #  Extract the output spectrum filename...
    outputFileName = args.outputFileName
    if args.verbose > 0:
        print 'outputFileName: ', outputFileName

    tab = ascii.read(inputFileName, format='csv')
    wave = tab['wave']  # Second column
    flux = tab['flux']  # Third column
    sp = S.ArraySpectrum(wave=wave, flux=flux, waveunits='angstroms', fluxunits='fnu')    

    sp.convert('flam')

    data = Table([sp.wave, sp.flux], names=['wave', 'flux'])
    ascii.write(data, outputFileName, format='csv')
    
    return 0

#--------------------------------------------------------------------------

if __name__ == "__main__":
    main()

#--------------------------------------------------------------------------

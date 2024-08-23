#!/usr/bin/env python
"""
    calc_abmag.py

    Calculate abmags for a spectrum (wavelength in Angstroms and flux in ergs/sec/cm**2/Angstrom 
    or in ergs/sec/cm**2/Hz) for a set of filter bandpass responses.

    The spectrum can be either in Synphot-style FITS format, or in CSV format.  If in CSV format, 
    one should designate the names for the wavelength column and the flux column.

    The filter bandpass responses should be in CSV format.  Examples can be found here: 
      https://github.com/DouglasLeeTucker/DECam_PGCM/tree/master/data/bandpasses


    Examples:

    calc_abmag.py --help

    calc_abmag.py --bandList g,r,i,z,Y --bandpassFile DES_STD_BANDPASSES_Y3A2_ugrizY.test.csv --spectrumFile SSSJ0006-5346_sum_bestfit.csv --colname_wave wave --colname_flux flux --flux_type Flam --verbose 1

    """


# Initial setup...
import numpy as np
import pandas as pd
import math
from scipy import interpolate
from astropy.io import fits
from astropy.table import Table
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

    parser.add_argument('--bandList', help='comma-separated list with no spaces', default='g,r,i,z,Y')
    parser.add_argument('--bandpassFile', help='name of the input plan file', default='DES_STD_BANDPASSES_Y3A2_ugrizY.test.csv')
    parser.add_argument('--spectrumFile', help='name of the input plan file (can be CSV file or a synphot-style FITS file')
    parser.add_argument('--colname_wave', help='name of the wavelength column (in case of a CSV spectrumFile)', default='wave')
    parser.add_argument('--colname_flux', help='name of the flux column (in case of a CSV spectrumFile)', default='flux')
    parser.add_argument('--flux_type', help='type of flux (Flam [ergs/sec/cm**2/Angstrom] or Fnu [ergs/sec/cm**2/Hz])? ', default='Flam')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    status = calc_abmag(args)

    return status

#--------------------------------------------------------------------------

def calc_abmag(args):

    #  Extract the bandList...
    bandList = args.bandList
    bandList = bandList.split(',')
    if args.verbose > 0:
        print 'bandList: ', bandList    
        
    #  Extract the name of the bandpassFile...
    bandpassFile = args.bandpassFile
    if os.path.isfile(bandpassFile)==False:
        print """bandpassFile %s does not exist...""" % (bandpassFile)
        print 'Returning with error code 1 now...'
        return 1
    if args.verbose > 0:
        print 'bandpassFile: ', bandpassFile

    #  Extract the name of the spectrum file...
    spectrumFile = args.spectrumFile
    if os.path.isfile(spectrumFile)==False:
        print """spectrumFile %s does not exist...""" % (spectrumFile)
        print 'Returning with error code 1 now...'
        return 1
    if args.verbose > 0:
        print 'spectrumFile: ', spectrumFile

    # Try to determine spectrumFile type (FITS file or CSV file)...
    spectrumType = 'Unknown'
    try:
        hdulist = fits.open(spectrumFile)
        hdulist.close()
        spectrumType = 'FITS'
    except IOError:
        if args.verbose > 2:
            print """spectrumFile %s is not a FITS file...""" % (spectrumFile)
        try:
            df_test = pd.read_csv(spectrumFile)
            spectrumType = 'CSV'
        except IOError:
            if args.verbose > 2:
                print """spectrumFile %s is not a CSV file...""" % (spectrumFile)

    # Read in spectrumFile and create a SciPy interpolated function of the spectrum...
    if spectrumType is 'FITS':
        flux,wave_lo,wave_hi = getSpectrumSynphot(spectrumFile, fluxFactor=1.0)
    elif spectrumType is 'CSV':
        flux,wave_lo,wave_hi = getSpectrumCSV(spectrumFile, fluxFactor=1.0)
    else:
        print """Spectrum file %s is of unknown type...""" % (spectrumFile)
        print 'Returning with error code 1 now...'
        return 1

    # Read the bandpassFile into a Pandas DataFrame...
    df_resp = pd.read_csv(bandpassFile, comment='#')

    # Check to make sure the spectrumFile covers at least the same wavelength range
    #  as the bandpassFile...
    if ( (wave_lo > df_resp['LAMBDA'].min()) or (wave_hi < df_resp['LAMBDA'].max()) ):
        print """WARNING:  %s does not cover the full wavelength range of %s""" % (spectrumFile, bandpassFile)
        print 'Returning with error code 1 now...'
        return 1

    # Create wavelength_array and flux_array...
    delta_wavelength = 1.0 # angstroms
    wavelength_array = np.arange(wave_lo, wave_hi, delta_wavelength)
    flux_array = flux(wavelength_array)

    
    # If needed, convert flux from flam to fnu...
    if args.flux_type == 'Fnu':
        fnu_array = flux_array
    elif args.flux_type == 'Flam':
        c_kms = 299792.5        # speed of light in km/s
        c_ms = 1000.*c_kms      # speed of light in m/s
        c_as = (1.000e10)*c_ms  # speed of light in Angstroms/sec
        fnu_array = flux_array * wavelength_array * wavelength_array / c_as
    else:
        print """Flux type %s is unknown...""" % (args.flux_type)
        print 'Returning with error code 1 now...'
        return 1


    # Print out header...
    outputLine = ''
    for band in bandList:
        outputLine =  """%s,%s""" % (outputLine, band)
    print outputLine[1:]

    outputLine = ''
    for band in bandList:

        response = interpolate.interp1d(df_resp['LAMBDA'], df_resp[band], 
                                        bounds_error=False, fill_value=0., 
                                        kind='linear')
        response_array = response(wavelength_array)

        try:
            abmag = calc_abmag_value(wavelength_array, response_array, fnu_array)
        except Exception:
            abmag = -9999.

        outputLine =  """%s,%.4f""" % (outputLine, abmag)

    print outputLine[1:]

    return 0


#--------------------------------------------------------------------------

# Calculate abmag using the wavelength version of the Fukugita et al. (1996) equation...

def calc_abmag_value(wavelength_array, response_array, fnu_array):

    # Calculate the abmag...
    numerator = np.sum(fnu_array * response_array / wavelength_array)
    denominator = np.sum(response_array / wavelength_array)
    abmag_value = -2.5*math.log10(numerator/denominator) - 48.60

    return abmag_value

#--------------------------------------------------------------------------

# Return a SciPy interpolation function of a Synphot-style FITS spectrum...
#  (Based on code from Keith Bechtol's synthesize_locus.py.)
# Unless otherwise noted, fluxes are assumed to be Flam and wavelengths  
#  are assumed to be in Angstroms...

def getSpectrumSynphot(synphotFileName, fluxFactor=1.0):

    try:

        hdulist = fits.open(synphotFileName)
        t = Table.read(hdulist[1])
        hdulist.close()

    except IOError:

        print """Could not read %s""" % synphotFileName
        sys.exit(1)


    wave = t['WAVELENGTH'].data.tolist()
    wave_lo = min(wave)
    wave_hi = max(wave)
    t['FLUX'] = fluxFactor*t['FLUX']
    flam = t['FLUX'].data.tolist()   
    flam = t['FLUX'].data.tolist()   
    data = {'wavelength': wave, 'flux': flam}

    f = interpolate.interp1d(data['wavelength'], data['flux'], 
                             bounds_error=True, 
                             kind='linear')

    return f,wave_lo,wave_hi

#--------------------------------------------------------------------------

# Return a SciPy interpolation function of a CSV-style spectrum...
#  (Based on code from Keith Bechtol's synthesize_locus.py.)
# Unless otherwise noted, fluxes are assumed to be Flam and wavelengths  
#  are assumed to be in Angstroms...

def getSpectrumCSV(csvFileName, colname_wave='wave', colname_flam='flux', fluxFactor=1.0):

    try:
        
        df = pd.read_csv(csvFileName)
        
    except IOError:

        print """Could not read %s""" % csvFileName
        sys.exit(1)


    columnNameList = df.columns.tolist()

    if colname_wave not in columnNameList:
        print """Column %s not in %s""" % (colname_wave, csvFileName)
        sys.exit(1)

    if colname_flam not in columnNameList:
        print """Column %s not in %s""" % (colname_wave, csvFileName)
        sys.exit(1)

    wave = df[colname_wave].tolist()
    wave_lo = min(wave)
    wave_hi = max(wave)
    df[colname_flam] = fluxFactor*df[colname_flam]
    flam = df[colname_flam].tolist()   
    data = {'wavelength': wave, 'flux': flam}

    f = interpolate.interp1d(data['wavelength'], data['flux'], 
                             bounds_error=True, 
                             kind='linear')

    return f,wave_lo,wave_hi

#--------------------------------------------------------------------------

# NOT CURRENTLY USED:

# Calculates and returns the Fitzpatrick 1999 reddening law
#  (for Rv=3.1) in inverse Angstroms...
# Based on code from Keith Bechtol's synthesize_locus.py code...

def getReddening_Fitz99():

    # Fitzpatrick 1999
    wavelength, a = zip(*[[2600, 6.591],
                          [2700, 6.265],
                          [4110, 4.315],
                          [4670, 3.806],
                          [5470, 3.055],
                          [6000, 2.688],
                          [12200, 0.829],
                          [26500, 0.265],
                          [1000000, 0.]])
    
    r = interpolate.interp1d(1. / np.array(wavelength), a, 
                             bounds_error=False, fill_value=0., kind=3) # NORMAL

    return r


#--------------------------------------------------------------------------

if __name__ == "__main__":
    main()

#--------------------------------------------------------------------------

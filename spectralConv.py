import pandas as pd
import numpy as np

def spectral_conv(R, S, x):
    '''
    returns the spectral convolution of spectral data, R with S, the spectral response function
    x are the wavelengths, which need to match both R and S (R and S should be same length and on same wavelength grid, x)
    see https://opg.optica.org/oe/fulltext.cfm?uri=oe-28-9-13801&id=431131 for details
    '''
    out = np.trapz(R * S, x=x) / np.trapz(S, x=x)
    return out

def SD_band_calc(RSR_doves, R, x):
    '''
    does spectral convolution for every superDoves band (although this is likely compatible with other sensors SRS too)
    RSR_doves is the spectral response data for superDoves (cols = Wavelength (nm), Coastal-Blue response, Blue etc..)
    R is the reflectance data to convolve for each band
    x is the wavelength grid for R and each spectral response function in RSR_doves
    '''
    R_SD = []
    for col in RSR_doves.columns[1:]: # don't use first col as this is the wavelengths column
        col_SR = RSR_doves[col].values
        R_SD.append(spectral_conv(R, col_SR, x))
    return np.array(R_SD)
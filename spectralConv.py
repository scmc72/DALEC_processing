import pandas as pd
import numpy as np
import dalecLoad

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

def SD_Rrs(RSR_doves, DALEC_sample, spect_wavelengths, x=None, doves_wavelengths=None, RHO=0.028, nsteps=601):
    '''
    does SD band calc for Lu, Lsky and Ed for a given DALEC sample, then converts this to Rrs using RHO
    returns a df with Lu, Lsky, Ed and Rrs
    '''
    Lu = dalecLoad.uniform_grid_spectra(DALEC_sample, spect_wavelengths, param='Lu', nsteps=nsteps)[:, 1]
    Lsky = dalecLoad.uniform_grid_spectra(DALEC_sample, spect_wavelengths, param='Lsky', nsteps=nsteps)[:, 1]
    Ed = dalecLoad.uniform_grid_spectra(DALEC_sample, spect_wavelengths, param='Ed', nsteps=nsteps)[:, 1]
    
    if x is None:
        x = RSR_doves['Wavelength (nm)'].values
    
    Lw = Lu - (RHO * Lsky)
    
    Lw_SD = SD_band_calc(RSR_doves, Lw, x)
    Lu_SD = SD_band_calc(RSR_doves, Lu, x)
    Lsky_SD = SD_band_calc(RSR_doves, Lsky, x)
    Ed_SD = SD_band_calc(RSR_doves, Ed, x)
    Rrs_SD = Lw_SD / Ed_SD
    
    if doves_wavelengths is None:
        doves_wavelengths = [443, 490, 531, 565, 610, 665, 705, 865]
    
    df_out = pd.DataFrame(data={'Wavelength': doves_wavelengths,
                           'Lu': Lu_SD, 
                           'Lw': Lw_SD,
                           'Lsky': Lsky_SD,
                           'Ed': Ed_SD,
                           'Rrs': Rrs_SD})
    
    return df_out


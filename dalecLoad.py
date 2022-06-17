import pandas as pd
import numpy as np

def load_DALEC_spect_wavelengths(filepath, header=15):
    """
    Loads spectral wavelength mappings from a calibrated DALEC logfile
    """
    spect_wavelengths = pd.read_csv(filepath,
                                    header=header,
                                    nrows=200,
                                    names=['Pixel_no', 'Ed', 'Lu', 'Lsky'],
                               )
    return spect_wavelengths # these are just the mappings of wavelength to pixel number

def load_DALEC_log(filepath, header=216, dropNA=True, longFormat=True, integerIndex=True):
    """
    loads DALEC log file (excluding spectral wavelength mappings)
    optionally returns log file in long format
    option to convert sample no. index to an integer, or to keep as a string (integerIndex)
    """
    DALEC_log = pd.read_csv(filepath,
                            header=header,
                            parse_dates=True,
                            dayfirst=True,
                            infer_datetime_format=True,
                            )
    if dropNA:
        DALEC_log.dropna(inplace=True, axis=0,)
        
    if longFormat:
        # convert to long format
        # need to test that these variable names always load in this way (leading space on Spec etc.)
        DALEC_log = pd.wide_to_long(DALEC_log, [' Spec['], i=['Sample #', ' Channel'], j='spectral_ind', suffix='\d+]')
        DALEC_log.reset_index(level=2, inplace=True) # remove spectral_ind as an index
        DALEC_log['spectral_ind'] = pd.to_numeric(DALEC_log['spectral_ind'].str[:-1]) # convert spectral_ind to numeric
        DALEC_log.rename(columns={' Spec[': 'Spectral Magnitude'}, inplace=True)
        DALEC_log = DALEC_log.astype({'Spectral Magnitude': 'float64'})
        if integerIndex:
            # change sample no. index to integer
            print('WARNING: some of my old code wont work with integerIndex=True - delete this line once this is sorted')
            idx = DALEC_log.index
            DALEC_log.index = DALEC_log.index.set_levels([idx.levels[0].astype(int), idx.levels[1]])
        # sort index
        DALEC_log.sort_index(inplace=True)

    return DALEC_log

from scipy import interpolate
# there are other interpolation methods too, but I think this is probably fine?

def uniform_grid_spectra(DALEC_sample, spect_wavelengths, param='Lu', nsteps=200, min_waveL=400, max_waveL=1000):
    """
    - takes spectrum from a single sample of a DALEC log file and converts to a uniform grid
    - grid is defined by nsteps, min_waveL and max_waveL
    - param gives which variable to grid: can choose between 'Lu', 'Lsky' and 'Ed' 
    """
    wavelength_grid = np.linspace(min_waveL, max_waveL, num=nsteps)
    
    y = DALEC_sample.loc[param]['Spectral Magnitude'].values
    x = spect_wavelengths[param].values
    interp = interpolate.interp1d(x, y)
    out = np.column_stack((wavelength_grid,
                         interp(wavelength_grid)))
    
    return out

def uniform_grid_spectra_mean(DALEC_log, spect_wavelengths, nsteps=601, min_waveL=400, max_waveL=1000):
    """
    - takes mean spectrum from an entire DALEC log file and converts to a uniform grid
    - grid is defined by nsteps, min_waveL and max_waveL
    - returns a pandas DF with Lu_mean, Lsky_mean and Ed_mean
    """
    # this is the bad code way to do it!
    Lu_tot = np.zeros((nsteps,))
    Lsky_tot = np.zeros((nsteps,))
    Ed_tot = np.zeros((nsteps,))
    for sample in DALEC_log.index.levels[0]:
        sample_i = DALEC_log.loc[sample, :]
        Lu_tot += uniform_grid_spectra(sample_i, spect_wavelengths, param='Lu', nsteps=nsteps)[:, 1]
        Lsky_tot += uniform_grid_spectra(sample_i, spect_wavelengths, param='Lsky', nsteps=nsteps)[:, 1]
        Ed_tot += uniform_grid_spectra(sample_i, spect_wavelengths, param='Ed', nsteps=nsteps)[:, 1]

    Lu_mean = Lu_tot/nsteps
    Lsky_mean = Lsky_tot/nsteps
    Ed_mean = Ed_tot/nsteps
    
#     # better way like this, which takes mean then does spectra regridding stuff
      # currently not working! need to work on this some more haha
#     DALEC_log.mean(axis=0, level='Sample #', numeric_only=True, inplace=True)
#     Lu_mean = dalecLoad.uniform_grid_spectra(DALEC_log, spect_wavelengths, param='Lu', nsteps=nsteps)
#     Lsky_mean = dalecLoad.uniform_grid_spectra(DALEC_log, spect_wavelengths, param='Lsky', nsteps=nsteps)
#     Ed_mean = dalecLoad.uniform_grid_spectra(DALEC_log, spect_wavelengths, param='Ed', nsteps=nsteps)
    
    wavelengths = uniform_grid_spectra(sample_i, spect_wavelengths, param='Lu', nsteps=nsteps)[:, 0]
    df_out = pd.DataFrame(data={'Wavelength': wavelengths,
                               'Lu_mean': Lu_mean, 
                               'Lsky_mean': Lsky_mean,
                               'Ed_mean': Ed_mean})
    return df_out

    
    
    
    
    
    
    
    
import pandas as pd
import numpy as np
import spectralConv
from scipy import interpolate
import os

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

def load_DALEC_log(filepath, header=216, dropNA=True, longFormat=True, integerIndex=True,
                   removeSaturated=True, parse_dates=[[' UTC Date', ' UTC Time']]):
    """
    loads DALEC log file (excluding spectral wavelength mappings)
    optionally returns log file in long format
    option to convert sample no. index to an integer, or to keep as a string (integerIndex)
    """
    # ideally specify dtype of all rows for efficiency and to prevent bad things - TODO!
    # need to specify str for lots of columns as these have some rows which contain stuff we need to remove
    DALEC_log = pd.read_csv(filepath,
                            header=header,
                            parse_dates=parse_dates,
                            dayfirst=True,
                            infer_datetime_format=True,
                            dtype={'Sample #': str,
                                   ' Lat': str, 
                                   ' Lon': str,
                                   ' Solar Azi': str,
                                   ' Solar Elev': str,
                                   ' Relaz': str,
                                   ' Heading': str,
                                   ' Pitch': str,
                                   ' Roll': str,
                                   ' Gearpos': str,
                                   ' Voltage': str,
                                   ' Temp': str,
                                   'Channel': str,
                                   ' Integration Time': str,
                                   ' Saturation Flag': str,
                                   ' Spec[21]': str,
                                   ' Spec[22]': str,
                                  },
                            )
    
    # any row with invalid UTC date can be removed
    DALEC_log.drop(DALEC_log[DALEC_log[' UTC Date_ UTC Time'].isna()].index, inplace = True)
    # this removes the duplicated headings
    DALEC_log.drop(DALEC_log[DALEC_log[' UTC Date_ UTC Time'] == ' UTC Date_ UTC Time'].index, inplace = True)
    
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
        # change saturation flag to int.
        DALEC_log[' Saturation Flag'] = DALEC_log[' Saturation Flag'].astype(int)
        # format column as datetimes
        DALEC_log[' UTC Date_ UTC Time'] = pd.to_datetime(DALEC_log[' UTC Date_ UTC Time'], dayfirst=True, infer_datetime_format=True)
        # remove saturated readings - this hasn't been tested on a df which isn't in long format!
        if removeSaturated:
            indSat = DALEC_log[DALEC_log[' Saturation Flag'] == 1].index.get_level_values(0)
            if list(indSat): # checks if the list is empty
                DALEC_log.drop(indSat, level=0, axis=0, inplace=True)
    DALEC_log.set_index(['spectral_ind', ' UTC Date_ UTC Time'], drop=True, append=True, inplace=True)
    
    # let's sort out these stupid names and remove leading spaces
    rename_cols_dict = dict([(col, col[1:]) for col in DALEC_log.columns[1:-1]])

    DALEC_log.rename(columns=rename_cols_dict, inplace=True)
    DALEC_log.index.rename(['Sample #', 'Channel', 'spectral_ind', 'Datetime'], inplace=True)
    
    # don't need sample no. as an index I don't think
    DALEC_log.index = DALEC_log.index.droplevel('Sample #')

    return DALEC_log

def load_DALEC_dir(DALEC_directory, file_names=None, **kwargs):
    '''
    finds all .dtf files in the specified directory dalec_dir and loads them using load_DALEC_log()
    results are combined into a single dataframe
    can also specify specific files in a list using file_names 
    '''
    if file_names is None: # if None, then load all DALEC transect (.dtf) files in the directory
        DALEC_files = []
        for file in os.listdir(DALEC_directory):
            if file.endswith(".dtf"):
                DALEC_files.append(os.path.join(DALEC_directory, file))
    else:
        DALEC_files = [DALEC_directory + file for file in file_names]
    
    print('loading ... ' + str(DALEC_files[0]))
    DALEC_df = load_DALEC_log(DALEC_files[0], **kwargs)
    
    for file in DALEC_files[1:]:
        print('loading ... ' + str(file))
        DALEC_df = pd.concat([DALEC_df, load_DALEC_log(file, **kwargs)])
            
    # probably smart to sort in case we get some weird stuff happenin' with file order etc.
    return DALEC_df.sort_values('Datetime')
    
    

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

def resampleMultiLog(dalecLog, freq='1D', level='Datetime', method='median'):
    '''
    - takes a dalecLog df (eg. logfile from load_DALEC_log())
    - resamples to the desired frequency (eg. daily '1D', or hourly '1h' etc.)
    - the aggregation method can be specified as a string, and defaults to median
    - any pandas aggregation function can be used (eg. 'mean', 'sum', 'std', 'min', 'max', 'describe')
    - most aggregation will only work on numeric data, and automatically removes non-numeric columns
    - (currently most cols are non-numeric, as that seems to be the easiest way to load em idk?)
    '''
    groupby = (dalecLog.groupby(['Channel', 'spectral_ind']
                                +[pd.Grouper(freq=freq, level='Datetime')]))
                           
    groupedMethod = getattr(groupby, method) # this basically allows you to call any method as a string
    return groupedMethod()
    
def uniform_grid_spectra_multi(DALEC_log, spect_wavelengths=None, RHO=0.028, nsteps=601, min_waveL=400, max_waveL=1000,
                               resample_to_SDs=True, col_end='_median'):
    '''
    takes a DALEC logfile (eg. from load_DALEC_log, or aggregated with resampleMultiLog)
    and performs regridding followed by calculation of Rrs 
    resample_to_SDs options allows for resampling to the superDoves wavebands
    NOTE this will regrid for all unique datetimes in the Datetime column and will be VERY SLOW
    if you have too many dates! hence, good idea to use resampleMultiLog() beforehand
    use col_end to add a suffix to the column names to indicate how they were previously calc'd
    eg. '_mean'
    '''
    
    if spect_wavelengths is None:
        spect_wavelengths = load_DALEC_spect_wavelengths('data/Jul-Aug/DALEC_72_73.dtf')
    
    df_out = None
    if resample_to_SDs:
        RSR_doves_file='non-DALEC-data/RSR-Superdove.csv'
        RSR_doves = pd.read_csv(RSR_doves_file)
        doves_wavelengths = [444., 492., 533., 566., 612., 666., 707., 866.]

    for date in DALEC_log.index.get_level_values('Datetime').unique():
        df = DALEC_log.loc[:, :, [date]]
        Lu = uniform_grid_spectra(df, spect_wavelengths, param='Lu', nsteps=nsteps)
        Lsky = uniform_grid_spectra(df, spect_wavelengths, param='Lsky', nsteps=nsteps)
        Ed = uniform_grid_spectra(df, spect_wavelengths, param='Ed', nsteps=nsteps)
        Rrs = (Lu[:, 1] - (RHO * Lsky[:, 1])) / Ed[:, 1]
        
        if resample_to_SDs:
            DALEC_SD = spectralConv.SD_band_calc(RSR_doves, Rrs,
                                                 RSR_doves['Wavelength (nm)'].values)

            df_tmp = pd.DataFrame(data=DALEC_SD, columns=['Rrs'+col_end])
            df_tmp['Date'] = np.full((8,), date)
            df_tmp['Wavelength'] = doves_wavelengths
            df_tmp.set_index(['Date', 'Wavelength'], inplace=True)

        else:
            df_tmp = pd.DataFrame(index=np.full((nsteps,), date),
                                  data={'Wavelength': Lu[:, 0],
                                        'Lu'+col_end: Lu[:, 1], 
                                        'Lsky'+col_end: Lsky[:, 1],
                                        'Ed'+col_end: Ed[:, 1],
                                        'Rrs'+col_end: Rrs})
            df_tmp.index.rename('Date', inplace=True)
            df_tmp.set_index('Wavelength', append=True, inplace=True)
            
        if df_out is None:
            df_out = df_tmp.copy()
        else:
            df_out = pd.concat([df_out, df_tmp])

    return df_out


def uniform_grid_spectra_mean(DALEC_log, spect_wavelengths, RHO=0.028, nsteps=601, min_waveL=400, max_waveL=1000):
    """
    - takes mean spectrum from an entire DALEC log file and converts to a uniform grid
    - grid is defined by nsteps, min_waveL and max_waveL
    - returns a pandas DF with Lu_mean, Lsky_mean and Ed_mean
    """
    print('for more flexibility and handling of multiple days its better to use resampleMultiLog() followed by uniform_grid_spectra_multi()')
    # this is a much more optimal way than previously! - 
    df = DALEC_log.copy() # not sure if neccesary but perhaps best to be on the safe side?
    # setting spectral_ind as an index might be useful for other stuff too?
    df = df.groupby(level=['Channel', 'spectral_ind']).mean(numeric_only=True)
    Lu_mean = uniform_grid_spectra(df, spect_wavelengths, param='Lu', nsteps=nsteps)
    Lsky_mean = uniform_grid_spectra(df, spect_wavelengths, param='Lsky', nsteps=nsteps)
    Ed_mean = uniform_grid_spectra(df, spect_wavelengths, param='Ed', nsteps=nsteps)
    
    Rrs_mean = (Lu_mean[:, 1] - (RHO * Lsky_mean[:, 1])) / Ed_mean[:, 1]
    
    df_out = pd.DataFrame(data={'Wavelength': Lu_mean[:, 0],
                               'Lu_mean': Lu_mean[:, 1], 
                               'Lsky_mean': Lsky_mean[:, 1],
                               'Ed_mean': Ed_mean[:, 1],
                               'Rrs_mean': Rrs_mean})
    return df_out

def uniform_grid_spectra_stats(DALEC_log, spect_wavelengths, RHO=0.028, nsteps=601, min_waveL=400, max_waveL=1000, 
                               percentiles=[.25, .5, .75],
                               fastGridding=True):
    """
    - finds summary stats from an entire DALEC log file and converts to a uniform grid
    - grid is defined by nsteps, min_waveL and max_waveL
    - returns a pandas DF with Lu_mean, Lsky_mean and Ed_mean
    """
    df = DALEC_log.copy() # not sure if neccesary but perhaps best to be on the safe side?
    
    # drop saturation flag to prevent this being included in the summary
    df.drop(labels='Saturation Flag', axis=1, inplace=True)
    df.set_index('spectral_ind', append=True, inplace=True)
    
    if fastGridding:
        print('WARNING: fastGridding enabled! - this will produce results much faster,'
              + ' but works by calculating Rrs before Lu, Lsky, and Ed have been interpolated'
              + ' to the same wavelength grid. Therefore, Rrs calculation may be inaccurate.')
        # previously needed to drop the Channel level, but now seems like not required... weird
        Lu = df.loc[:, 'Lu', :]['Spectral Magnitude']#.droplevel(' Channel')
        Lsky = df.loc[:, 'Lsky', :]['Spectral Magnitude']#.droplevel(' Channel')
        Ed = df.loc[:, 'Ed', :]['Spectral Magnitude']#.droplevel(' Channel')
    
        df_Rrs = (Lu - (RHO * Lsky)) / Ed
        df_Rrs = df_Rrs.groupby(level=['spectral_ind']).describe(percentiles=percentiles)
        wavelength_grid = np.linspace(min_waveL, max_waveL, num=nsteps)
    
        y = df_Rrs.values
        x = spect_wavelengths['Lu'].values

        interp = interpolate.interp1d(x, y, axis=0)
        Rrs_summary = np.column_stack((wavelength_grid,
                                   interp(wavelength_grid)))
        colnames = ['Wavelength'] + list(df_Rrs.columns) # get column names
        df_out = pd.DataFrame(data=Rrs_summary, columns=colnames)

    else:
        print('WARNING fast gridding disabled. Interpolation of Lu, Lsky and Ed may be VERY SLOW \
               for datasets with large numbers of samples. Try enabling fastGridding for much faster,\
               but potentially less accurate results... \n \n \
               NOTE: "slowGridding" not yet developed. \n \n \
               RETURNING None')

        df_out = None
    
    #df = df.groupby(level=[' Channel', 'spectral_ind']).describe(percentiles=percentiles)
    

    #Lu_summary = uniform_grid_spectra(df, spect_wavelengths, param='Lu', nsteps=nsteps)
    # this should have shape (nsteps, 9), where columns are:
    # wavelength, count, mean, std, min, 25%, 50%, 75%, max
    # note that cols will be different if length of percentiles arg changes
    #Lsky_summary = uniform_grid_spectra(df, spect_wavelengths, param='Lsky', nsteps=nsteps)
    #Ed_summary = uniform_grid_spectra(df, spect_wavelengths, param='Ed', nsteps=nsteps)
    
    # performing this calculation with the mean makes sense, but I'm not sure it actually makes sense
    # for the SD, median etc... NEED TO FIX THIS
    
    # options: calc Rrs before regridding and calc'ing percentiles etc. (will be fast)
    # regrid entire dataset and calc Rrs for everything, then do summary of this directly (will be SLOWWWW)
    
    #Rrs_summary = (Lu_summary[:, 1:] - (RHO * Lsky_summary[:, 1:])) / Ed_summary[:, 1:]
    #Rrs_summary = uniform_grid_spectra(df_Rrs, spect_wavelengths, param='Rrs', nsteps=nsteps)
    
    #colnames = ['wavelength'] + [x[1] for x in list(df.columns)] # get column names
    #data = np.concatenate([Lu_summary[:, 0].reshape(nsteps, 1), Rrs_summary], axis=1)
#    df_out = pd.DataFrame(data=data, columns=colnames)
    #df_out = pd.DataFrame(data=Rrs_summary, columns=colnames)
    return df_out

def uniform_grid_spectra_Rrs(DALEC_sample, spect_wavelengths, RHO=0.028, nsteps=601, min_waveL=400, max_waveL=1000):
    '''
    takes a single sample from a DALEC log and does spectrum gridding followed by basic Rrs calculation
    '''
    Lu = uniform_grid_spectra(DALEC_sample, spect_wavelengths, param='Lu', nsteps=nsteps)[:, 1]
    Lsky = uniform_grid_spectra(DALEC_sample, spect_wavelengths, param='Lsky', nsteps=nsteps)[:, 1]
    Ed = uniform_grid_spectra(DALEC_sample, spect_wavelengths, param='Ed', nsteps=nsteps)[:, 1]
    Rrs = (Lu - (RHO * Lsky)) / Ed
    
    wavelengths = uniform_grid_spectra(DALEC_sample, spect_wavelengths, param='Lu', nsteps=nsteps)[:, 0]
    
    df_out = pd.DataFrame(data={'Wavelength': wavelengths,
                               'Lu': Lu, 
                               'Lsky': Lsky,
                               'Ed': Ed,
                               'Rrs': Rrs})
    return df_out


def multiLogLoad(filepath, 
                 sep=['DALEC (SN:0005)'],
                 header=216, 
                 dropNA=True,longFormat=True, 
                 integerIndex=True,
                 removeSaturated=True):
    """
    - loads multiple logs which are all contained in a single logfile (eg. that was generated using serial logging of DALEC
    - depending on how the logfile was generated, adjusting 'sep' might allow for different situations...
    - should probably think a bit more about this for a continuous logging application
    """
    # ideally specify dtype of all rows for efficiency and to prevent bad things - TODO!
    # need to specify str for lots of columns as these have some rows which contain stuff we need to remove
    
    print('probably dont need to use this function anymore as load_DALEC_log() should support loading of these file types now!')
    
    df = pd.read_csv(filepath,
                     header=header,
                     parse_dates=True,
                     dayfirst=True,
                     infer_datetime_format=True,
                     dtype={'Sample #': str,
                                   ' Lat': str, 
                                   ' Lon': str,
                                   ' Solar Azi': str,
                                   ' Solar Elev': str,
                                   ' Relaz': str,
                                   ' Heading': str,
                                   ' Pitch': str,
                                   ' Roll': str,
                                   ' Gearpos': str,
                                   ' Voltage': str,
                                   ' Temp': str,
                                   'Channel': str,
                                   ' Integration Time': str,
                                   ' Saturation Flag': str,
                                   ' Spec[21]': str,
                                   ' Spec[22]': str,
                                  },
                     )
    
    groups = df['Sample #'].isin(sep).cumsum()
    names = ['Log ' + str(i) for i in range(len(set(groups)))] 
    tables = {name: g[1].iloc[1:] for g,name in zip(df.groupby(groups), names)} 
    # because we've used 'DALEC (SN:0005)' as the way to seperate logfiles, we need to remove the first lines before this
    # if a different sep is used, then this line may no longer work... so perhaps need to make this more robust!
    tables.pop('Log 0')
    
    for name, table in tables.items():
        # any row with invalid UTC date can be removed
        table.drop(table[table[' UTC Date'].isna()].index, inplace = True)
        # this removes the duplicated headings
        table.drop(table[table[' UTC Date'] == 'UTC Date'].index, inplace = True)
        if dropNA:
            table.dropna(inplace=True, axis=0,)
        if longFormat:
            # convert to long format
            # need to test that these variable names always load in this way (leading space on Spec etc.)
            table = pd.wide_to_long(table, [' Spec['], i=['Sample #', ' Channel'], j='spectral_ind', suffix='\d+]')
            table.reset_index(level=2, inplace=True) # remove spectral_ind as an index
            table['spectral_ind'] = pd.to_numeric(table['spectral_ind'].str[:-1]) # convert spectral_ind to numeric
            table.rename(columns={' Spec[': 'Spectral Magnitude'}, inplace=True)
            table = table.astype({'Spectral Magnitude': 'float64'})
        if integerIndex:
            # change sample no. index to integer
            print('WARNING: some of my old code wont work with integerIndex=True - delete this line once this is sorted')
            idx = table.index
            table.index = table.index.set_levels([idx.levels[0].astype(int), idx.levels[1]])
        # sort index
        table.sort_index(inplace=True)
        # change saturation flag to int.
        table[' Saturation Flag'] = table[' Saturation Flag'].astype(int)
        # format column as datetimes
        table[' UTC Date'] = pd.to_datetime(table[' UTC Date'], dayfirst=True, infer_datetime_format=True)
        # remove saturated readings - this hasn't been tested on a df which isn't in long format!
        if removeSaturated:
            indSat = table[table[' Saturation Flag'] == 1].index.get_level_values(0)
            if list(indSat): # checks if the list is empty
                table.drop(indSat, level=0, axis=0, inplace=True)
        
        tables[name] = table
    return tables
    

    
    
    
    
    
    
    
    
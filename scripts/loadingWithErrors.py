import pandas as pd

# this is now pretty good!
# I can load multiple .dtf files from a directory, including those which span multiple days
# I can resample very easily and generate some good plots with the data


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
        
    print(DALEC_log)
    
    for i in range(88574, 89348):
        
        if len(DALEC_log[DALEC_log['Sample #']==str(i)]) != 3:
            print('error at sample # = ' + str(i))
            print(DALEC_log[DALEC_log['Sample #']==str(i)])
    DALEC_log.drop_duplicates(subset=['Sample #', ' Channel'], keep='last', inplace=True)
        
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

                        
data_path = 'C:/Users/daa5/Project/DALEC_processing/data/Jul-Aug/LOG_0090.dtf'
                      
loglog = load_DALEC_log(data_path, longFormat=True)

print(loglog)



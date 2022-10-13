import dalecLoad as dl
import pandas as pd
import SD_NC_loading as sdl

# this is now pretty good!
# I can load multiple .dtf files from a directory, including those which span multiple days
# I can resample very easily and generate some good plots with the data

print('loading...')

# multiLog = dl.load_DALEC_log('data/Jul-Aug/DALEC_72_73.dtf')

# result = dl.resampleMultiLog(multiLog, method='median')

# df_daily_gridded = dl.uniform_grid_spectra_multi(result, resample_to_SDs=False)
# print(df_daily_gridded)

log_71_73 = dl.load_DALEC_dir('data/Jul-Aug/')#, file_names=['DALEC_71.dtf', 'DALEC_72_73.dtf'])
print('log_71_73:')
print(log_71_73)

result = dl.resampleMultiLog(log_71_73, method='median', freq='10min')
df_daily_gridded = dl.uniform_grid_spectra_multi(result, resample_to_SDs=False)
print('daily median vals gridded:')
print(df_daily_gridded)

#sdl.multiDaySpectraPlot(df_daily_gridded,
                        # SD_col_slice=slice(0, 0), # slice returns empty as we don't have any SD data yet
                        # ylim=[0, df_daily_gridded['Rrs_median'].max(axis=0)*1.1],
                        # figsize=(12, 12))
                        
# df_daily_SD = dl.uniform_grid_spectra_multi(result, resample_to_SDs=True)
# sdl.multiDaySpectraPlot(df_daily_SD,
                        # SD_col_slice=slice(0, 0), 
                        # figsize=(12, 12))
                        
                        
sdl.plot_algorithm_from_DF(df_daily_gridded, col_names=['Rrs_median'])

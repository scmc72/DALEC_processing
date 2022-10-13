# idea is to run this script followed by Dalin's code to get DALEC data into a nice format ready for processing 

import os
import sys
lib_path = os.path.abspath(os.path.join(os.path.abspath(''), '../functions/'))
print(lib_path)
sys.path.append(lib_path)
import dalecLoad as dl
import pandas as pd
import SD_NC_loading as sdl
import time

import numpy as np
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt

outputFileName = 'Smoothin_and_groovin.csv' # perhaps would be good to have some naming convention?
savWindowLen = 21 # 21 nm is recommeneded by Jiang et al., 2020, ISPRS.
                  # https://doi.org/10.1016/j.isprsjprs.2020.05.003
savPolyOrder = 4 # Dalin doesn't specify in the paper what order, but I think he uses R and 4 appears to be the deafult  

data_path = os.path.abspath(os.path.join(os.path.abspath(''), '../data/'))

log = dl.load_DALEC_dir(os.path.join(data_path, 'Jul-Aug'))#,#, # choose data folder here
                        #file_names=['LOG_0090.dtf']) # remove file_names to load all files in directory

# can adjust resampling frequency and method here:
result = dl.resampleMultiLog(log, method='median', freq='10min')

df_daily_gridded = dl.uniform_grid_spectra_multi(result, resample_to_SDs=False)

df_daily_gridded.drop(labels=['Lu_median', 'Lsky_median', 'Ed_median'], axis=1, inplace=True)
df1 = df_daily_gridded.unstack().add_prefix('Rrs') # get into wide format
df1 = df1['RrsRrs_median'] # remove abstract col header Rrs_median
df1.columns = [col[:-2] for col in df1.columns] # remove '.0' from end of col names

# print(df1)

# arr = df1.to_numpy()
# print('before smoothin:')
# print(arr)
# print('array shape = ' + str(arr.shape))

# arrSmooth = savgol_filter(arr, savWindowLen, savPolyOrder, axis=1)
# print('after smoothin:')
# print(arrSmooth)

# APPLY SAVITZKY GOLAY FILTER  
dfSmooth = df1.apply(lambda x: savgol_filter(x, savWindowLen, savPolyOrder),
                     axis=1, result_type='broadcast')
                     
# arrSmooth = dfSmooth.to_numpy()

# plt.plot(arr[0, :])
# plt.plot(arrSmooth[0, :])
# plt.plot(arr[12, :])
# plt.plot(arrSmooth[12, :])
# plt.show()
          
          

Rcode_path = os.path.abspath(os.path.join(os.path.abspath(''), '../R-code/'))

outputPath = os.path.join(Rcode_path, 'data/' + outputFileName)

if os.path.exists(outputPath):
    print('file already exists, adding current time to avoid overwriting...')
    outputPath = outputPath[:-4] + str(time.time()) + '.csv'
    
dfSmooth.to_csv(outputPath)
print('saved: ' + outputPath)

# save as .csv ready for Dalin's delta correction code
           
           

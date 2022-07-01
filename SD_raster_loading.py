import numpy as np
import rasterio
import rasterio.sample 
from rasterio.plot import reshape_as_raster, reshape_as_image
from pyproj import Transformer
import pandas as pd
import matplotlib.pyplot as plt

def getSpectraFromSDSR(rasterFile, xx, yy, crs_SD='epsg:32630', crs_coords= 'WGS 84'):
    '''
    - takes a filename to a superdoves surface reflectance tif file and xx, yy coords as inputs
    xx, yy can be lists, arrays or scalars
    - returns a pandas df with columns for each coord and the SD wavelengths
    - default: take input coords in WGS 84 format (N, W: eg. 56.147209, -3.923337) and then index the tif using epsg:32630
    but this can be changed where needed (eg. if looking at data outside UK then epsg:32630 might not be appropriate)
    CRS must be specified in the pyproj string format 
    '''
    transformer = Transformer.from_crs(crs_coords, crs_SD)
    xxT, yyT = transformer.transform(xx, yy)
    with rasterio.open(rasterFile) as dataset:
        try:
            spect = rasterio.sample.sample_gen(dataset, zip(xxT, yyT))
        except TypeError:
            spect = rasterio.sample.sample_gen(dataset, zip([xxT], [yyT]))
        wavelengths = pd.DataFrame(data={'Band Name':dataset.descriptions,
                                         'Wavelength':[443, 490, 531, 565, 610, 665, 705, 865]})
        data = np.array(list(spect)).T/(2**16)
    df = pd.DataFrame(data=data)
    df = wavelengths.join(df)
    return df

def getSpectraFromSDSR_grid(rasterFile, x, y, shape=(3, 3), crs_SD='epsg:32630', crs_coords= 'WGS 84'):
    '''
    generates grid with shape = shape around the coord (x, y) of interest (default crs for coord WGS84)
    returns spectra for each coord in pandas df
    '''
    transformer = Transformer.from_crs(crs_coords, crs_SD)
    xxT, yyT = transformer.transform(x, y)
    
    with rasterio.open(rasterFile) as dataset:
        row, col = dataset.index(xxT, yyT)
        # generate x and y coords for grid with shape=(shape[0], shape[1])
        rows = np.linspace(row - shape[0]//2,
                        row + shape[0]//2 - (1 - shape[0]%2),
                        shape[0],
                        dtype=int)
        cols = np.linspace(col - shape[1]//2,
                        col + shape[1]//2 - (1 - shape[1]%2),
                        shape[1],
                        dtype=int)
      
        inds = (np.tile(rows, len(cols)), np.repeat(cols, len(rows)))
        spect = dataset.read()[np.index_exp[:] + inds]
        wavelengths = pd.DataFrame(data={'Band Name':dataset.descriptions,
                                         'Wavelength':[443, 490, 531, 565, 610, 665, 705, 865]})
        
        data = np.array(list(spect))/(2**16)
        
    df = pd.DataFrame(data=data)
    df = wavelengths.join(df)
    return df

def plotSDRaster(rasterFile, ax=None, overDrive=1.0, plotShow=False):
    '''
    quick function to get a plot of the specified superdoves raster file
    overDrive is used to add brightness to the image by multiplying the RGB values (which have been scaled to 0-1) by this value
    overDrive > 1.0 will introduce some clipping to the image displayed
    '''
    with rasterio.open(rasterFile) as dataset:
        data = dataset.read()
    image = reshape_as_image(data)
    # this slice selects bands 2, 4, and 6 (blue, green, red), then flip to get RGB in correct order
    image = np.flip(image[:, :, 1:6:2], 2).astype(np.float64)
    # excellent use of for loop for normalisation (0.0->1.0)
    for i in range(3):
        image[:, :, i] = (image[:, :, i] - image[:, :, i].min()) / (image[:, :, i].max() - image[:, :, i].min())
        
    if ax is None:
        ax = plt.gca()
    ax.imshow(image * overDrive)
    if plotShow:
        plt.show()
        
def plotSDRasterSpectraGrid(rasterFile, x, y, ax=None, shape=(3, 3), crs_SD='epsg:32630', crs_coords= 'WGS 84'):
    '''
    plots spectra from superdoves raster file for a grid of pixels, shape=shape, at a given x, y point
    (default CRS is WGS84, but can change using crs_coords)
    '''
    df = getSpectraFromSDSR_grid(rasterFile, x, y, shape=(3, 3), crs_SD='epsg:32630', crs_coords= 'WGS 84')
    
    if ax is None:
        ax = plt.gca()
    
    for col in list(df.columns.values)[2:]:
        ax.plot(df['Wavelength'],
                 df[col],
                 color='green',
                 label='SuperDoves SR',
                 marker='o',
                 alpha=0.2)
        
    handles, labels = ax.get_legend_handles_labels()
    newLabels, newHandles = [], []
    for handle, label in zip(handles, labels):
        if label not in newLabels:
            newLabels.append(label)
            newHandles.append(handle)
    ax.legend(newHandles, newLabels)
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('$R_{rs}$ $(sr^{-1}$)')    
    ax.grid()
    plt.show()
            
            
            
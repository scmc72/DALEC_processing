import numpy as np
from scipy.signal import savgol_filter
x = np.array([[2, 2, 5, 2, 1, 0, 1, 4, 9], [2, 2, 5, 2, 1, 0, 1, 4, 9]]).T

print(x.shape)
a = savgol_filter(x, 5, 2, axis=0)

print(a)
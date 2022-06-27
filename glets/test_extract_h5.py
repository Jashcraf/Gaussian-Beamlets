import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

f = h5.File('/Users/jashcraft/Desktop/gbd-data/M3.h5', 'r')

keys = f.keys()
print(keys)
meas = f['measurement0']
print(meas.keys())
dat = meas['reserve_interferogram']['frame4']['data']
# print(dat.keys())



plt.figure()
plt.imshow(dat) # vmin=-1e-1,vmax=5e-1
plt.colorbar()
plt.show()
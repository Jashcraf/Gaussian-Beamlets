import numpy as np
import glets.utils as ut

fn = '/Users/jashcraft/Desktop/gbd-data/Hubble_RTM.mat'
coordsfn = '/Users/jashcraft/Desktop/gbd-data/pupilcoords.txt'

coords = np.loadtxt(coordsfn,delimiter=',')
x = coords[:,0]
y = coords[:,1]
print(coords.shape)

data = ut.ReadMatlabMat(fn)
dim = data.shape[-1]
dim = int(np.sqrt(dim))
x = np.linspace(-1,1,dim)
x,y = np.meshgrid(x,x)

ut.fourbyfour(data,x,y,size=1)
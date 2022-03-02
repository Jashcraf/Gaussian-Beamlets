#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 10:32:33 2021

@author: jashcraft
"""
import numpy as np
import numpy
import matplotlib.pyplot as plt
import matplotlib as mpl
from skimage.restoration import unwrap_phase

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.size'] = 18

def grab_center(array,cut):
    
    d1,d2 = array.shape
    return array[d1-cut:d1+cut,d2-cut:d2+cut]

def amp_phase_plot(array,logintensity=True):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    
    plt.set_cmap('plasma')
    fig,ax = plt.subplots(ncols=2,figsize=[14,7])
    if logintensity == True:
        cmapi = ax[0].imshow(np.abs(array),norm=LogNorm(),origin='lower')
    else:
        cmapi = ax[0].imshow(np.abs(array),origin='lower')
    
    
    ax[0].set_title('Irradiance')
    fig.colorbar(cmapi,ax=ax[0])
        
    plt.set_cmap('coolwarm')
    cmapp = ax[1].imshow(np.angle(array),origin='lower')
    ax[1].set_title('Phase')
    fig.colorbar(cmapp,ax=ax[1])
    
    plt.show()
    
def angularspectrum(array,pixelscale,wavelength):
    
    spectrum = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(array)))
    spectrum_phase = unwrap_phase(np.angle(spectrum))
    spectrum_amp   = np.abs(spectrum)
    
    spectrum_angle_x,spectrum_angle_y = np.gradient(spectrum_phase,pixelscale)
    
    return spectrum,wavelength*spectrum_angle_x,wavelength*spectrum_angle_y

def hexagonal_grid(radius,spacing=1.):
    
    cosv = np.cos(np.pi/6)
    sinv = np.sin(np.pi/6)
    
    nsteps = int(radius/(spacing*cosv))
    i = np.arange(-nsteps-1, nsteps+2)
    j = np.arange(-nsteps-1, nsteps+2)
    
    vi,vj = np.meshgrid(i,j)
    
    x = (vi + sinv*vj)*spacing
    y = cosv*vj*spacing
    
    r2 = x**2 + y**2
    select = r2 < (radius**2)
    
    return x[select], y[select]

def fourbyfour(array,size,coords=None):
    
    import matplotlib.tri as tri
    
    # What I think they should be
    # amax = 1e-3
    # amin = -amax
    # bmax = 100
    # bmin = -bmax
    # cmax = 1/bmax #1/bmax
    # cmin = -1/bmax #-1/bmax
    # dmin = -1.1
    # dmax = 1.1
    
    amax = 1#1e-2 #None
    amin = 0#-1e-2 #None
    bmax = 1e-6#254e-3 # 5.52e-5 #None
    bmin = 0 # -5.52e-5 #None
    cmax = -1/255e-3 #1/5.52e-3 # None
    cmin = -1/254e-3 #Non e #min(array[2,0,:])
    dmax = 1e-0#1e-4 # None
    dmin = 0 # None
    
    # amax = None
    # amin = None
    # bmax = None
    # bmin = None
    # cmax = None
    # cmin = None 
    # dmax = None
    # dmin = None
    
    Axx = array[0,0,:]
    Axy = array[0,1,:]
    Ayx = array[1,0,:]
    Ayy = array[1,1,:]
    
    Bxx = array[0,2,:]
    Bxy = array[0,3,:]
    Byx = array[1,2,:]
    Byy = array[1,3,:]
    
    Cxx = array[2,0,:]
    Cxy = array[2,1,:]
    Cyx = array[3,0,:]
    Cyy = array[3,1,:]
    
    Dxx = array[2,2,:]
    Dxy = array[2,3,:]
    Dyx = array[3,2,:]
    Dyy = array[3,3,:]
    
    if coords is not None:
        x = coords[0,:]
        y = coords[1,:]
        print(x.shape)
    abcd = [Axx,Axy,Ayx,Ayy,Bxx,Bxy,Byx,Byy,Cxx,Cxy,Cyx,Cyy,Dxx,Dxy,Dyx,Dyy]
    titl = ['Axx','Axy','Ayx','Ayy','Bxx','Bxy','Byx','Byy','Cxx','Cxy','Cyx','Cyy','Dxx','Dxy','Dyx','Dyy']
    
    # for i in range(16):
        
    #     # plt.figure()
    #     # plt.imshow(np.reshape(abcd[i],[50,50]),vmin=0,vmax=1e-5)
    #     # plt.title(titl[i])
    #     # plt.colorbar()
    #     # plt.show()
        
    #     plt.figure()
    #     if coords is not None:
    #         plt.scatter(x,y,c=abcd[i])
    #     else:
    #         x = np.linspace(-size/2,size/2,int(np.sqrt(array[0,0,:].size)))
    #         x,y = np.meshgrid(x,x)
    #         x = np.ravel(x)
    #         y = np.ravel(y)
    #         plt.scatter(x,y,c=abcd[i])
    #     plt.title(titl[i])
    #     plt.colorbar()
    #     plt.show()
        
    
    fig,ax = plt.subplots(ncols=4,nrows=4,figsize=[20,16])
    
    plt.suptitle('Ray Transfer Matrix')
    
    pca = ax[0,0].scatter(x,y,c=Axx,vmin=amin,vmax=amax)
    ax[0,0].axis('off')
    ax[0,1].scatter(x,y,c=Axy,vmin=amin,vmax=amax)
    ax[0,1].axis('off')
    ax[1,0].scatter(x,y,c=Ayx,vmin=amin,vmax=amax)
    ax[1,0].axis('off')
    ax[1,1].scatter(x,y,c=Ayy,vmin=amin,vmax=amax)
    ax[1,1].axis('off')
    
    pcb = ax[0,2].scatter(x,y,c=Bxx,vmin=bmin,vmax=bmax)
    ax[0,2].axis('off')
    ax[0,3].scatter(x,y,c=Bxy,vmin=bmin,vmax=bmax)
    ax[0,3].axis('off')
    ax[1,2].scatter(x,y,c=Byx,vmin=bmin,vmax=bmax)
    ax[1,2].axis('off')
    ax[1,3].scatter(x,y,c=Byy,vmin=bmin,vmax=bmax)
    ax[1,3].axis('off')
    
    pcc = ax[2,0].scatter(x,y,c=Cxx,vmin=cmin,vmax=cmax)
    ax[2,0].axis('off')
    ax[2,1].scatter(x,y,c=Cxy,vmin=cmin,vmax=cmax)
    ax[2,1].axis('off')
    ax[3,0].scatter(x,y,c=Cyx,vmin=cmin,vmax=cmax)
    ax[3,0].axis('off')
    ax[3,1].scatter(x,y,c=Cyy,vmin=cmin,vmax=cmax)
    ax[3,1].axis('off')
    
    pcd = ax[2,2].scatter(x,y,c=Dxx,vmin=dmin,vmax=dmax)
    ax[2,2].axis('off')
    ax[2,3].scatter(x,y,c=Dxy,vmin=dmin,vmax=dmax)
    ax[2,3].axis('off')
    ax[3,2].scatter(x,y,c=Dyx,vmin=dmin,vmax=dmax)
    ax[3,2].axis('off')
    ax[3,3].scatter(x,y,c=Dyy,vmin=dmin,vmax=dmax)
    ax[3,3].axis('off')
    
    fig.colorbar(pca, ax=ax[0:2,0:2], shrink=0.6, location='right')
    fig.colorbar(pcb, ax=ax[0:2,2:4], shrink=0.6, location='right')
    fig.colorbar(pcc, ax=ax[2:4,0:2], shrink=0.6, location='right')
    fig.colorbar(pcd, ax=ax[2:4,2:4], shrink=0.6, location='right')
    
    plt.show()
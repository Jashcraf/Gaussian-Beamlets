#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 08:17:58 2021

@author: jashcraft
"""

import numpy as np
from scipy.interpolate import RectBivariateSpline
import poppy
import numexpr as ne
import time
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from scipy.ndimage import zoom
import utils as ju

def makerays(size,numrays,circle=True):
    
    # Define lists of XY coordinate pairs for square grid
    x = np.linspace(-size/2,size/2,numrays)
    y = np.linspace(-size/2,size/2,numrays)
    x,y = np.meshgrid(x,y)
    X = np.ravel(x)
    Y = np.ravel(y)
    r = np.sqrt(X**2 + Y**2)
    X = X[r <= size/2]
    Y = Y[r <= size/2]
    
    return np.array([X,Y,0*X,0*Y])

def eval_gausfield(rays,sys,wavelength,detsize,npix,amps=None):
    """
    

    Parameters
    ----------
    rays : ndarray
        4 x nrays array of base rays (wow).
    sys : ndarray
        4 x 4 x nrays ray transfer matrix list.
    wavelength : float
        wavelength of gaussians in meters.
    detsize : float
        dimension of detector in meters.
    npix : int
        number of pixels along dimension detsize.
    amps : ndarray
        list of gaussian amplitudes.

    Returns
    -------
    ndarray.

    """
    
    # Prepare propagation bits
    wo = 10*wavelength
    zr = np.pi*wo**2/wavelength
    Q = np.array([[1/(1j*zr),0],[0,1/(1j*zr)]])
    
    u = np.linspace(-detsize/2,detsize/2,npix)
    u,v = np.meshgrid(u,u)
    u = np.ravel(u)
    v = np.ravel(v)
    Dphase = np.zeros([len(u),len(rays[0,:])],dtype='complex')
    
    if amps == None:
        amps = 1 + 0*rays[0,:]
    
    for i in np.arange(0,rays.shape[1]):
        
        # Step 1 - Propagate Rays
        rayp = np.matmul(sys[:,:,i],rays[:,i])
        
        A = sys[0:2,0:2,i]
        B = sys[0:2,2:4,i]
        C = sys[2:4,0:2,i]
        D = sys[2:4,2:4,i]
        
        # Step 2 - Propagate Complex Beam Parameter
        Qp_n = (C + np.matmul(D,Q))
        Qp_d = np.linalg.inv(A+np.matmul(B,Q))
        Qp   = np.matmul(Qp_n,Qp_d)
        
        if A[0,0] == 0:
            orig_matrx = np.zeros([2,2])
        else:
            orig_matrx = np.linalg.inv(np.linalg.inv(Q) + np.matmul(np.linalg.inv(A),B))
            
        cros_matrx = np.linalg.inv(np.matmul(A,np.linalg.inv(Q))+B)
        
        lo = sys[0,2,i]
        
        # Step 2 - Evaluate Phasor
        uo = u-rays[0,i]
        vo = v-rays[1,i]

        up = u-rayp[0]
        vp = v-rayp[1]

        # And this is where all the math comes in (Lin et al)
        guoy_phase = -1j*np.arctan(lo/np.real(Qp[0,0]))
        tran_phase = (-1j*(np.pi/wavelength))*(Qp[0,0]*up**2 + (Qp[1,0] + Qp[0,1])*up*vp + Qp[1,1]*vp**2)
        long_phase = -1j*(2.0*np.pi/wavelength)*lo

        orig_phase = (-1j*(np.pi/wavelength))*(orig_matrx[0,0]*uo**2 + (orig_matrx[1,0] + orig_matrx[0,1])*uo*vo + orig_matrx[1,1]*vo**2)

        cros_phase = (-1j*(2*np.pi/wavelength))*( cros_matrx[0,0]*uo*up + (cros_matrx[1,0] + cros_matrx[0,1])*uo*vp + cros_matrx[1,1]*vo*vp )
        Dphase[:,i] = tran_phase+long_phase+guoy_phase+orig_phase+cros_phase
    
    t1 = time.perf_counter()
    Efield = np.sum(amps*ne.evaluate('exp(Dphase)'),axis=1)
    t2 = time.perf_counter()
    print('Time to compute exponential ',t2-t1,' s')
    
    return Efield

def IdentityMat(nrays):
    
    box = np.zeros([4,4,nrays])
    
    # construct identity matrix
    box[0,0,:] = 1.
    box[1,1,:] = 1.
    box[2,2,:] = 1.
    box[3,3,:] = 1.
    
    return box
    

def AnamorphicLens(eflx,efly,nrays):
    
    analens = IdentityMat(nrays)
    analens[2,0,:] = -1/eflx
    analens[3,1,:] = -1/efly
    
    return analens

def ThinLens(efl,nrays):
    
    return AnamorphicLens(efl,efl,nrays)

def FreeSpace(distance,nrays):
    
    dmat = IdentityMat(nrays)
    dmat[0,2,:] = distance
    dmat[1,3,:] = distance
    
    return dmat

def OPDSurfaceError(nrays,pixscal):
    
    # note: need to find where rays hit on the surface 
    # RectBivariateSpline()
    npix = 512
    x = np.linspace(-1,1,npix)
    x,y = np.meshgrid(x,x)
    zoomby = int(np.sqrt(nrays))/512
    
    opdmat = IdentityMat(nrays)
    
    #opd = kolphase(512)*1e-8
    #opd = zernike_nm(5,3,512)*1e-9
    opd = (x**2 - y**3)*1e-7
    
    # it appears that I need to add an imperceptible difference in order to
    # not mask opd when I mask opdmask? i'm not sure why
    opdmask = opd + 1e-30
    
    
    opdmask[x**2 + y**2 >= 1] = 0
    
    # do some zero-padding
    osbox = np.zeros([4096,4096])
    osbox[2048-256:2048+256,2048-256:2048+256] = opdmask
    
    
    psfopd = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(np.exp(1j*2*np.pi*osbox/633e-9))))
    cut = 32
    psfopd = psfopd[2048-cut:2048+cut,2048-cut:2048+cut]
    
    plt.figure(figsize=[7,7])
    plt.imshow(opd)
    plt.colorbar()
    plt.title('OPD Applied in meters')
    plt.show()
    
    ju.amp_phase_plot(psfopd)
    
    gradx,grady = np.gradient(zoom(opd,zoomby),pixscal)
    
    opdmat[2,0,:] = np.ravel(gradx)
    opdmat[3,1,:] = np.ravel(grady)
    
    return opdmat
    
def kolphase(s):
    #phase=kolphaseClass(s)
    # Returns a random Kolmogorov phase screen of dimension s x s computed 
    # from the FT of random complex numbers with appropriate amplitudes. 
    # Screens are computed on a grid of size 2s, with a s x s piece cut out. 
    # This helps overcome the problem with this techniqe of under-representing 
    # tilt. Tandom tilts are also explicity added to give a reasonable 
    # approximation of the overall Kolmogorav structure function on all scales.
    ph=np.zeros([s,s,2]) #initialize phase variable
    [x,y]=np.meshgrid(range(-s,s),range(-s,s))
    r=np.sqrt(np.multiply(x,x)+np.multiply(y,y)) #make a radial ordinate
    
    f1=np.random.randn(2*s,2*s) #make two sets of gaussian random numbers
    f2=np.random.randn(2*s,2*s)
    
    f=f1+1j*f2 #turn the two random numbers into a complex random number
    ps=np.power(r,-11/6) #Kolomogorov power spectrum amplitude
    ps[s][s]=0; #setting the centr of the power spectrum (DC term) to zero
    
    scale=15.2 #this number makes the tilts work out 
    
    xt=np.random.randn(2)*(scale/s) #random amplitudes for tip and tilt to 
    yt=np.random.randn(2)*(scale/s) #mitigate power spectral rolloff at low frequencies 
    
    #This has Kolmogorov phases in real and imaginary parts
    sc= np.fft.fft2(np.fft.fftshift(np.multiply(np.multiply(np.abs(f),ps),np.exp(1j*np.angle(f))))) 
    
    #seperating the real and imaginary parts
    ph[:,:,0]=np.real(sc[0:s,0:s])+xt[0]*x[0:s,0:s]+yt[0]*y[0:s,0:s]
    ph[:,:,1]=np.real(sc[0:s,0:s])+xt[1]*x[0:s,0:s]+yt[1]*y[0:s,0:s]
    
    #make the mean phase zero
    ph[:,:,0]=ph[:,:,0]-np.mean(ph[:,:,0])
    ph[:,:,1]=ph[:,:,1]-np.mean(ph[:,:,1])
    
    #pick just one of the two phase screens for present purposes scale to 
    #something that will give reasonable ansewrs on a 1024 grid with 256 pupil
    ph=ph[:,:,1]*3
    
    return ph

    
def OPD2GridSag(opdarray,pixscal):
    
    dim = opdarray.shape[0]
    unitflag = 3 # in meters
    xdec = 0
    ydec = 0
    
    
    row1 = [dim,dim,pixscal,pixscal,unitflag,xdec,ydec]
    
    box = np.zeros([opdarray.size,5])
    box[:,0] = np.ravel(opdarray)*1e-7
    
    np.savetxt('kolphase.DAT',box)
    
def matmultlist(mat1,mat2):
    
    # this works for 4x4xn dimension arrays, where it multiplies the matrices element-wise by the final axis
    
    box = np.zeros(mat1.shape)
    
    for mind in range(mat1[0,0,:].size):
        
        box[:,:,mind] = np.matmul(mat1[:,:,mind],mat2[:,:,mind])
        
    return box

def zernikeRadialFunc(n, m, r):
    """
    Fucntion to calculate the Zernike radial function

    Parameters:
        n (int): Zernike radial order
        m (int): Zernike azimuthal order
        r (ndarray): 2-d array of radii from the centre the array

    Returns:
        ndarray: The Zernike radial function
    """
    import numpy
    R = numpy.zeros(r.shape)
    for i in range(0, int((n - m) / 2) + 1):

        R += numpy.array(r**(n - 2 * i) * (((-1)**(i)) *
                         numpy.math.factorial(n - i)) /
                         (numpy.math.factorial(i) *
                          numpy.math.factorial(0.5 * (n + m) - i) *
                          numpy.math.factorial(0.5 * (n - m) - i)),
                         dtype='float')
    return R

def zernike_nm(n, m, N):
    """
     Creates the Zernike polynomial with radial index, n, and azimuthal index, m.

     Args:
        n (int): The radial order of the zernike mode
        m (int): The azimuthal order of the zernike mode
        N (int): The diameter of the zernike more in pixels
     Returns:
        ndarray: The Zernike mode
     """
    import numpy 
    coords = (numpy.arange(N) - N / 2. + 0.5) / (N / 2.)
    X, Y = numpy.meshgrid(coords, coords)
    R = numpy.sqrt(X**2 + Y**2)
    theta = numpy.arctan2(Y, X)

    if m==0:
        Z = numpy.sqrt(n+1)*zernikeRadialFunc(n, 0, R)
    else:
        if m > 0: # j is even
            Z = numpy.sqrt(2*(n+1)) * zernikeRadialFunc(n, m, R) * numpy.cos(m*theta)
        else:   #i is odd
            m = abs(m)
            Z = numpy.sqrt(2*(n+1)) * zernikeRadialFunc(n, m, R) * numpy.sin(m * theta)

    # clip
    Z = Z*numpy.less_equal(R, 1.0)

    return Z#*circle(N/2., N)


def proprays(rays,system):
    
    box = np.zeros(rays.shape)
    
    for rind in range(rays[0,:].size):
        
        box[:,rind] = np.matmul(system[:,:,rind],rays[:,rind])
        
    return box


a = time.perf_counter()

npix = 256
rays = makerays(25.4e-3,26,circle=True)
sys = matmultlist(OPDSurfaceError(rays.shape[1],25.4e-3/256),ThinLens(1,rays.shape[1]))
sys = matmultlist(FreeSpace(1,rays.shape[1]),sys)
field = np.reshape(eval_gausfield(rays,sys,633e-9,50e-5,npix),[npix,npix])

import utils as ju
ju.fourbyfour(sys,25.4e-3)


import matplotlib.pyplot as plt

ju.amp_phase_plot(field)


        
    
    
    
    
    
    
    
    
    
    
    
    


        
    
    # Step 3 - Compute Field
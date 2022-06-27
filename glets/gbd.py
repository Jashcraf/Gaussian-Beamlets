#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 12:57:22 2022

@author: jashcraft
"""

import numpy as np
import time
import multiprocessing
import os
os.environ['NUMEXPR_MAX_THREADS'] = '96'
os.environ['NUMEXPR_NUM_THREADS'] = '8'
import numexpr as ne
from numba import jit, set_num_threads
print('numexpr sees numcores = ',ne.ncores)
print('numexpr sees numthreads = ',ne.nthreads)
print('numexpr sees max threads = ',ne.MAX_THREADS)

numbathreads = 8
print('numba threads = ',numbathreads)
set_num_threads(numbathreads)

def Matmulvec(x2,y2,M,x1,y1):

    return (x2*M[0,0] + y2*M[1,0])*x1 + (x2*M[0,1] + y2*M[1,1])*y1

def Inv(M):
    a = M[0,0]
    b = M[0,1]
    c = M[1,0]
    d = M[1,1]

    return 1/(a*d - b*c)*np.array([[d,-b],
                                   [-c,a]])

def Det(M):
    a = M[0,0]
    b = M[0,1]
    c = M[1,0]
    d = M[1,1]

    return (a*d-b*c)

def Mult(M1,M2):

    a1 = M1[0,0]
    b1 = M1[0,1]
    c1 = M1[1,0]
    d1 = M1[1,1]

    a2 = M2[0,0]
    b2 = M2[0,1]
    c2 = M2[1,0]
    d2 = M2[1,1]

    A = a1*a2 + b1*c2
    B = a1*b2 + b1*d2
    C = c1*a2 + d1*c2
    D = c1*b2 + d1*d2

    mat = np.array([[A,B],[C,D]])

    return mat

def ComputeGouyPhase(Q):

    eigvals = np.linalg.eigvals(Q)
    q1,q2 = eigvals[0],eigvals[1]

    gouy = .5*(np.arctan(np.real(q1)/np.imag(q1)) + np.arctan(np.real(q2)/np.imag(q2)))

    return gouy

    

def EvalGausslets(x,y,rays,Qinv,k,sys,npix,use_numexpr=True,raysout=None,lo=0,use_numba=False):
    """Evaluate a coherent sum of decentered astigmatic gaussian beamlets

    Parameters
    ----------
    x : numpy.ndarray
        detector coordinate in the x direction
    y : numpy.ndarray
        detector coordinate in the y direction
    rays : 4 x 1 x nrays numpy.ndarray
        center rays = [x,y,u,v] of the gaussian at the FIRST surface of the optical system
    Qinv : 2x2 numpy.ndarray
        inverse of the complex curvature matrix
    k : scalar
        wavenumber = 2*pi/wavelength
    sys : 4 x 4 x nrays numpy.ndarray
        ray transfer matrices composed of 2x2 submatrices [A,B;C,D]
    npix : scalar
        number of resolution elements along one dimension of a square detector
    use_numexpr : bool, optional
        Use numexpr instead of numpy to evaluate the exponential, by default True
    raysout : 4 x 1 x nrays numpy.ndarray, optional
        center rays = [x,y,u,v] of the gaussian at the LAST surface of the optical system, 
        by default None - and is computed with the supplied ABCD
    lo : scalar, optional
        Eikonal along the optical axis, by default 0

    Returns
    -------
    numpy.ndarray
        Complex amplitude of the coherent sum of gaussian beams
    """

    # Get number of rays to trace
    nrays = rays.shape[1]
    print("Number of Gausslets = ",nrays)

    # Define empty box to put complex amplitudes in
    field = np.zeros([int(npix*npix)],dtype='complex128')

    # Consider Multiprocessing w/ Numba if it's too much, pretty slow right now
    for i in range(nrays):

        # Need to make sure that these are actually aligned! check against pupilcoords
        ray = rays[:,i]
        ABCD = sys[:,:,i]

        # This is a bad conditional - find an alternative
        if type(raysout) == np.ndarray:

            rout = raysout[:,i]

        else:

            rout = None

        if type(lo) == np.ndarray:

            loin = lo[i]

        else:
            loin = 0

        field += EvalGaussian(x,y,ray,Qinv,k,ABCD,use_numexpr=True,rayout=rout,lo=loin)

    print('finished tracing ',i,' Gausslets')

    field = np.reshape(field,[npix,npix])

    return field

def EvalGaussian(x,y,ray,Qinv,k,ABCD,use_numexpr=True,rayout=None,lo=None):
    """Evaluate complex amplitude of a decentered astigmatic gaussian

    Parameters
    ----------
    x : numpy.ndarray
        detector coordinate in the x direction
    y : numpy.ndarray
        detector coordinate in the y direction
    ray : 4x1 numpy.ndarray
        center ray = [x,y,u,v] of the gaussian at the FIRST surface of the optical system
    Qinv : 2x2 numpy.ndarray
        inverse of the complex curvature matrix
    k : scalar
        wavenumber = 2*pi/wavelength
    ABCD : 4x4 numpy.ndarray
        ray transfer matrix composed of 2x2 submatrices [A,B;C,D]
    use_numexpr : bool, optional
        Use numexpr instead of numpy to evaluate the exponential, by default True
    rayout : numpy.ndarray, optional
        center ray = [x,y,u,v] of the gaussian at the LAST surface of the optical system, 
        by default None - and is computed with the supplied ABCD
    lo : scalar, optional
        Eikonal along the optical axis of a ray, by default None

    Returns
    -------
    numpy.ndarray
        Complex amplitude of a DAGB at every point on the detector
    """

    phasor,amplitude = PropagateGaussianPhase(x,y,ray,Qinv,k,ABCD,rayout,lo)

    if use_numexpr == True:

        field = amplitude*ne.evaluate('exp(phasor)')

    else:

        field = amplitude*np.exp(phasor)

    return field

def PropagateGaussianPhase(x,y,ray,Qinv,k,ABCD,rayout,lo):
    """Compute phase of propagated decentered astigmatic gaussian beam (DAGB)

    Parameters
    ----------
    x : numpy.ndarray
        detector coordinate in the x direction
    y : numpy.ndarray
        detector coordinate in the y direction
    ray : 4x1 numpy.ndarray
        center ray = [x,y,u,v] of the gaussian at the FIRST surface of the optical system
    Qinv : 2x2 numpy.ndarray
        inverse of the complex curvature matrix
    k : scalar
        wavenumber = 2*pi/wavelength
    ABCD : 4x4 numpy.ndarray
        ray transfer matrix composed of 2x2 submatrices [A,B;C,D]
    rayout : 4x1 numpy.ndarray
        center ray = [x,y,u,v] of the gaussian at the LAST surface of the optical system
    lo : scalar
        Eikonal along the optical axis of the ray

    Returns
    -------
    numpy.ndarray
        phase of the DAGB evaluated at every point on the detector
    """

    # Grab Decenter Parameters
    xd = ray[0]
    yd = ray[1]

    if type(rayout) == np.ndarray:
        ray_prop = rayout
    else:
        ray_prop = ABCD @ ray

    # Grab submatrices of optical system
    A = ABCD[0:2,0:2]
    B = ABCD[0:2,2:4]
    C = ABCD[2:4,0:2]
    D = ABCD[2:4,2:4]

    # shift by beamlet position
    x = x-ray_prop[0]
    y = y-ray_prop[1]

    # Ravel detector coordinates for faster computing
    x = np.ravel(x)
    y = np.ravel(y)

    # compute the information from the matrices in Cai and Lin Eq 9
    Qpinv = (C + D @ Qinv) @ np.linalg.inv(A + B @ Qinv)
    Q = np.linalg.inv(Qinv)
    Qp = np.linalg.inv(Qpinv)

    # Try your own hacked inverse
    # Qpinv = (C + D @ Qinv) @ Inv(A + B @ Qinv)
#     Qpinv = Mult((C + Mult(D , Qinv)) , Inv(A + Mult(B , Qinv)))
#     Q = Inv(Qinv)

    # Phase from gaussian beam
    gaus = GaussianPhase(x,y,Qpinv,k)

    # Phase from decenter
    demat = np.linalg.inv(Q + np.linalg.inv(A) @ B)

    # Try your own hacked inverse
    # demat = Inv(Q + Inv(A) @ B)
#     demat = Inv(Q + Mult(Inv(A) , B))

    # Phase from cross term
    crmat = np.linalg.inv(A @ Q + B)

    # Try your own hacked inverse
#     crmat = Inv(Mult(A , Q) + B)


    # phase from decenter parameter
    # This is part of the aberration - but without this it looks less comatic
    gaus += -1j*k/2 * Matmulvec(xd,yd,demat,xd,yd) # (xd**2 * demat[0,0] + (demat[1,0] + demat[0,1])*xd*yd + demat[1,1]*yd**2)

    # phase from coupling of decenter and position
    # Without this it just looks gaussian
    gaus += 1j*k * Matmulvec(xd,yd,crmat,x,y) #(xd*crmat[0,0]*x + xd*crmat[1,0]*y + yd*crmat[0,1]*x + yd*crmat[1,1]*y)

    # apply eikonal along optical axis
    lo = 57.6
    gaus += 1j*k*lo

    # apply gouy phase shift
    gaus += 1j*ComputeGouyPhase(Qp)

    # amplitude scaling
    # This parameter complicates the PSF - see Issue #40
    # For now setting it to unity
    amp = 1/np.sqrt((np.linalg.det(A + B @ Qinv)))
#     amp = 1/np.sqrt((Det(A + Mult(B , Qinv))))

    return gaus,amp

def GaussianPhase(x1,y1,Qinv,k):
    """Determine phase of a gaussian beam

    Parameters
    ----------
    x1 : numpy.ndarray
        detector coordinate in the x direction
    y1 : numpy.ndarray
        detector coordinate in the y direction
    Qinv : 2x2 numpy.ndarray
        inverse of the omplex curvature matrix
    k : scalar
        wavenumber = 2*pi/wavelength

    Returns
    -------
    numpy.ndarray
        phase of an astigmatic gaussian beam on a detector plane
    """
    
    return -1j*k/2 * Matmulvec(x1,y1,Qinv,x1,y1)#(x1**2 * Qinv[0,0] + (Qinv[1,0] + Qinv[0,1])*x1*y1 + Qinv[1,1]*y1**2)

def MakeEvenGrid(size,wo,OF):

    numrays = np.int(np.round(size*OF/(2*wo)))

    # Define lists of XY coordinate pairs for square grid
    x = np.linspace(-size/2,size/2,numrays)
    y = np.linspace(-size/2,size/2,numrays)
    x,y = np.meshgrid(x,y)
    X = np.ravel(x)
    Y = np.ravel(y)

    return X,Y

def MakeFibonacciGrid(size,wo,OF):

    numrays = np.int(np.round(np.pi*((size/2.0)*OF/(wo))*9.0)) # This is arbitrary scaling
    print('numbeamlets = ',numrays)
        
    
    c = np.array([0,0]) # XY offset from a spiral
    R = (size/2)*np.sqrt(np.linspace(1/2,numrays-1/2,numrays))/np.sqrt(numrays-1/2)
    T = 4/(1+np.sqrt(5))*np.pi*np.linspace(1,numrays,numrays)
    X = c[0] +R*np.cos(T)
    Y = c[1] +R*np.sin(T)

    return X,Y

def ComputeRaysFromOF(OF,size,wo,samplescheme='even'):
    
    if samplescheme == 'fibbonacci':
        print('Fibonacci OF Not solved yet, accuracy may vary')
        X,Y = MakeFibonacciGrid(size,wo,OF)

    elif samplescheme == 'even':

        X,Y = MakeEvenGrid(size,wo,OF)
        
    return np.array([X,Y,0*X,0*Y])
        
        

def eval_gausfield(rays,sys,wavelength,detsize,npix,amps=None,use_numba=False,use_numexpr=False,return_rays=False):
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
    wo = 9200*wavelength
    zr = np.pi*wo**2/wavelength
    Q = np.array([[1/(1j*zr),0],[0,1/(1j*zr)]])
#     size = 2*np.abs(rays[0,0])
#     print('pupil D = ',size)
#     print('OF = ',2*wo*np.sqrt(rays.shape[1])/size)
    
    u = np.linspace(-detsize/2,detsize/2,npix)
    u,v = np.meshgrid(u,u)
    u = np.ravel(u)
    v = np.ravel(v)
    Dphase = np.zeros([len(u),len(rays[0,:])],dtype='complex128')
    
    if amps == None:
        amps = 1 + 0*rays[0,:]
    
    t00 = time.perf_counter()
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
        guoy_phase = -1j*np.pi/2 #np.arctan(lo/np.real(Qp[0,0]))
        tran_phase = (-1j*(np.pi/wavelength))*(Qp[0,0]*up**2 + (Qp[1,0] + Qp[0,1])*up*vp + Qp[1,1]*vp**2)
        long_phase = -1j*(2.0*np.pi/wavelength)*lo

        orig_phase = (-1j*(np.pi/wavelength))*(orig_matrx[0,0]*uo**2 + (orig_matrx[1,0] + orig_matrx[0,1])*uo*vo + orig_matrx[1,1]*vo**2)

        cros_phase = (-1j*(2*np.pi/wavelength))*( cros_matrx[0,0]*uo*up + (cros_matrx[1,0]*up*vo + cros_matrx[0,1]*uo*vp) + cros_matrx[1,1]*vo*vp )
        Dphase[:,i] = (tran_phase+long_phase+guoy_phase+orig_phase+cros_phase)
    
    t01 = time.perf_counter()
    print('time to assemble phase array = ',t01-t00)
    
    if use_numexpr == True:
        print('profiling numexpr')
        t1 = time.perf_counter()
        Efield = np.sum(ne.evaluate('exp(Dphase)'),axis=1) # amps*
        t2 = time.perf_counter()
    else:
        print('profiling numpy')
        t1 = time.perf_counter()
        Efield = np.sum(np.exp(Dphase),axis=1)
        t2 = time.perf_counter()
        
    print('Time to compute exponential ',t2-t1,' s')
    
    return Efield
    
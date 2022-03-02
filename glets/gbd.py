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
os.environ['NUMEXPR_NUM_THREADS'] = '96'
import numexpr as ne
from numba import jit, set_num_threads
print('numexpr sees numcores = ',ne.ncores)
print('numexpr sees numthreads = ',ne.nthreads)
print('numexpr sees max threads = ',ne.MAX_THREADS)

numbathreads = 16
print('numba threads = ',numbathreads)
set_num_threads(numbathreads)

def eval_gausfield(rays,sys,wavelength,detsize,npix,amps=None,use_numba=False,use_numexpr=False):
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
    wo = 25*wavelength
    zr = np.pi*wo**2/wavelength
    Q = np.array([[1/(1j*zr),0],[0,1/(1j*zr)]])
    
    u = np.linspace(-detsize/2,detsize/2,npix)
    u,v = np.meshgrid(u,u)
    u = np.ravel(u)
    v = np.ravel(v)
    Dphase = np.zeros([len(u),len(rays[0,:])],dtype='complex')
    
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

        cros_phase = (-1j*(2*np.pi/wavelength))*( cros_matrx[0,0]*uo*up + (cros_matrx[1,0] + cros_matrx[0,1])*uo*vp + cros_matrx[1,1]*vo*vp )
        Dphase[:,i] = (tran_phase+long_phase+guoy_phase+orig_phase+cros_phase)
    
    t01 = time.perf_counter()
    print('time to assemble phase array = ',t01-t00)
    
    if use_numexpr == True:
        print('profiling numexpr')
        t1 = time.perf_counter()
        Efield = np.sum(ne.evaluate('exp(Dphase)'),axis=1) # amps*
        t2 = time.perf_counter()
    elif use_numba == True:
        print('profiling numba')
        t1 = time.perf_counter()
        Efield = numba_eval_gausfield(Dphase)
        t2 = time.perf_counter()
    else:
        print('profiling numpy')
        t1 = time.perf_counter()
        Efield = np.sum(np.exp(Dphase),axis=1)
        t2 = time.perf_counter()
        
    print('Time to compute exponential ',t2-t1,' s')
    
    return Efield

@jit(nopython=True, parallel=True)
def numba_eval_gausfield(dphase):
    return np.sum(np.exp(dphase),axis=1)
    

def target_eval_gausfield(ray,sys,wavelength,detsize,npix,amps=None):
    
    # Prepare propagation bits
    wo = 10*wavelength
    zr = np.pi*wo**2/wavelength
    Q = np.array([[1/(1j*zr),0],[0,1/(1j*zr)]])
    
    u = np.linspace(-detsize/2,detsize/2,npix)
    u,v = np.meshgrid(u,u)
    u = np.ravel(u)
    v = np.ravel(v)
    Dphase = np.zeros([len(u)],dtype='complex')
    
    if amps == None:
        amps = 1 + 0*ray[0]
    
        
    # Step 1 - Propagate Rays
    rayp = np.matmul(sys,ray)
    
    A = sys[0:2,0:2]
    B = sys[0:2,2:4]
    C = sys[2:4,0:2]
    D = sys[2:4,2:4]
    
    # Step 2 - Propagate Complex Beam Parameter
    Qp_n = (C + np.matmul(D,Q))
    Qp_d = np.linalg.inv(A+np.matmul(B,Q))
    Qp   = np.matmul(Qp_n,Qp_d)
    
    if A[0,0] == 0:
        orig_matrx = np.zeros([2,2])
    else:
        orig_matrx = np.linalg.inv(np.linalg.inv(Q) + np.matmul(np.linalg.inv(A),B))
        
    cros_matrx = np.linalg.inv(np.matmul(A,np.linalg.inv(Q))+B)
    
    lo = sys[0,2]
    
    # Step 2 - Evaluate Phasor
    uo = u-ray[0]
    vo = v-ray[1]

    up = u-rayp[0]
    vp = v-rayp[1]

    # And this is where all the math comes in (Lin et al)
    guoy_phase = -1j*np.arctan(lo/np.real(Qp[0,0]))
    tran_phase = (-1j*(np.pi/wavelength))*(Qp[0,0]*up**2 + (Qp[1,0] + Qp[0,1])*up*vp + Qp[1,1]*vp**2)
    long_phase = -1j*(2.0*np.pi/wavelength)*lo

    orig_phase = (-1j*(np.pi/wavelength))*(orig_matrx[0,0]*uo**2 + (orig_matrx[1,0] + orig_matrx[0,1])*uo*vo + orig_matrx[1,1]*vo**2)

    cros_phase = (-1j*(2*np.pi/wavelength))*( cros_matrx[0,0]*uo*up + (cros_matrx[1,0] + cros_matrx[0,1])*uo*vp + cros_matrx[1,1]*vo*vp )
    Dphase[:,i] = tran_phase+long_phase+guoy_phase+orig_phase+cros_phase
    
    Efield = np.sum(np.exp(Dphase),axis=1)
    

if __name__ == "__main__":
    
    p1 = multiprocessing.target_eval_gausfield()
    
    
"""
Engine for doing gaussian beamlet propagation
"""
import numpy as np
import numexpr as ne 

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
    
    Efield = np.sum(amps*ne.evaluate('exp(Dphase)'),axis=1)
    
    return Efield

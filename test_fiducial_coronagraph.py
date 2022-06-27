# Beamlet scripts
import glets.raytrace as rt
import glets.gbd as gf
import glets.utils as ut

# Utilities
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.ndimage import zoom

# poppy
import astropy.units as u
import poppy

# Simple profiling
import time

from scipy.ndimage import shift
from matplotlib.colors import LogNorm

def ComputeIntensity(array):
    return np.abs(array*np.conj(array))

def array2poppyhdul(array,pixscal):
    
    from astropy.io import fits
    phdu = fits.PrimaryHDU(array)
    hdul = fits.HDUList([phdu])
    hdul[0].header['PIXELSCL'] = pixscal
    
    return hdul

def compute_mean_std(data):
    return np.mean(data),np.std(data)

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

def propagate_NIRCam(wfarray):
    
    """
    The IWA and OWA using an EFL = 132.812m the sizes are
    IWA = 0.01544m
    OWA = 0.04633m
    """
    
    osys = poppy.FresnelOpticalSystem(npix=256, beam_ratio=2)
    #osys.add_optic(poppy.ScalarTransmission())
    
    
    # The Annular Occulting Mask
    #     osys.add_optic(poppy.CircularAperture(radius=2.5e-5*7*u.m))
    # had to arbitrarily scale - I think we missed a factor of 10 somewhere?
    osys.add_optic(poppy.InverseTransmission(poppy.CircularAperture(radius=0.0001544*u.m)))
    osys.add_optic(poppy.CircularAperture(radius=0.0004633*u.m))
    
    # The Unit magnification relay 10" lens - efl is arbitrary since it's 1-1
    osys.add_optic(poppy.QuadraticLens(254e-3*u.m),distance=254e-3*u.m)
    
    # The Lyot Stop
    osys.add_optic(poppy.CircularAperture(radius=25.4e-3*0.9*u.m),distance=254e-3*u.m)
    
    osys.add_optic(poppy.QuadraticLens(254e-3*u.m),distance=254e-3*u.m)
    osys.add_optic(poppy.ScalarTransmission(),distance=254e-3*u.m)
    psf = osys.calc_psf(wavelength=2.2e-6*u.m,inwave=wfarray)
    
    return psf

def FiducialCoronagraphTest(error=None):

    # error corresponds to phase error amount
    

    # Set up Osys parameters for Hybrid and Fresnel System
    nrays_across = 101
    scale = 1/2
    wlen = 2.1e-6
    npix = 256
    size = 6.5
    efl = 132.812
    detector_size = 5.36356153808263e-06*npix#0.004633 #5e-4
    # poppy.conf.use_numexpr = True
    # poppy.accel_math._USE_NUMEXPR

    phase = kolphase(npix)*1e-8

    print('begin fresnel setup')
    wf = poppy.FresnelWavefront(size/2*u.m,wavelength=wlen,npix=npix,oversample=4)
    circ = poppy.CircularAperture(radius=size/2*u.m)
    pwfe = poppy.ArrayOpticalElement(transmission=np.ones(phase.shape),opd=phase,pixelscale=wf.pixelscale)
    amount = error # 3*2.1e-6
    # pwfe = poppy.ZernikeWFE(coefficients=[0,0,0,0+amount,0,0],radius=2*size/2*u.m)
    thin = poppy.QuadraticLens(efl*u.m)

    # Grab applied WFE
    opd_to_pass = ut.grab_center(pwfe.get_opd(wf),npix/2)

    # plt.figure(figsize=[10,5])
    # pwfe.display(what='both',opd_vmax=50e-7)
    # plt.show()

    print('propagating fresnel')

    wf *= pwfe
    wf *= circ
    wf *= thin

    z = efl*u.m # F/10 lens
    wf.propagate_fresnel(z)

    plt.figure(figsize=(12,5))
    wf.display('both',colorbar=True,imagecrop=detector_size,scale='log')
    plt.suptitle("After propagating {}".format(z), fontsize=18)
    plt.show()

    print('Begin GBD setup')

    # Create square equally-spaced grid of rays
    # rays = rt.makerays(size, nrays_across,circle=False)
    rays = gf.ComputeRaysFromOF(1.4,size,20000*wlen)
    print(rays.shape)
    nrays_total = int((rays.shape[1]))
    nrays_across = int(np.sqrt(nrays_total))

    # Give rays some wavefront error
    array_gbd = zoom(opd_to_pass,nrays_across/npix)
    testrays = rt.ArbitraryWFE(nrays_across, size,rays,scale,array=array_gbd) # array assumes units of meters

    # Set up optical system
    lens = rt.ThinLens(efl, nrays_total)
    dist = rt.FreeSpace(efl, nrays_total)
    rayt = rt.matmultlist(lens, dist)
    rayt = rt.matmultlist(dist,rayt)

    # Select rays from circular aperture
    raysincirc = testrays[:,np.sqrt(testrays[0,:]**2 + testrays[1,:]**2) <= size/2]
    rtmincirc  = rayt[:,:,np.sqrt(testrays[0,:]**2 + testrays[1,:]**2) <= size/2]

    print('propagate GBD')

    # evaluate coherent sum of gaussian beamlets
    efie_numerical = gf.eval_gausfield(raysincirc, rtmincirc, wlen, detector_size, npix,return_rays=True)
    efie_numerical = np.reshape(efie_numerical,[npix,npix])

    ut.amp_phase_plot(efie_numerical,logintensity=True)

    # grab PSF center
    efie_fresnel = ut.grab_center(wf.wavefront,npix/2)
    vmin = None
    vmax = None
    norm = LogNorm()

    fnlpsf = ComputeIntensity(efie_fresnel)
    fnlpsf /= np.sum(fnlpsf)
    gbdpsf = shift(ComputeIntensity(efie_numerical),[0.5,0.5])
    gbdpsf /= np.sum(gbdpsf)
    difpsf = fnlpsf-gbdpsf

    pixscale = 5.36356153808263e-06

    efie_numerical /= np.sum(efie_numerical)
    efie_numerical /= np.max(efie_numerical)
    efie_fresnel /= np.sum(efie_fresnel)
    efie_fresnel /= np.max(efie_fresnel)

    wfg = poppy.FresnelWavefront(beam_radius=5e-4/2*u.m,npix=256,oversample=2,wavelength=2.1e-6*u.m)
    wff = poppy.FresnelWavefront(beam_radius=5e-4/2*u.m,npix=256,oversample=2,wavelength=2.1e-6*u.m)
    wfg.wavefront = zoom(efie_numerical,2)
    wff.wavefront = zoom(efie_fresnel,2)

    psfg = propagate_NIRCam(wfg)
    psff = propagate_NIRCam(wff)
    psfd = ComputePSFDifference(psfg,psff,psfg[0].header['PIXELSCL'])

    # Compute Coronagraphic Profiles
    gbd_pro = poppy.radial_profile(psfg)
    fnl_pro = poppy.radial_profile(psff)
    dif_pro = gbd_pro[1]-fnl_pro[1]

    plt.figure(figsize=[22,15])
    plt.suptitle('Coronagraphic Focal Plane Comparison, RMS = {:.2E}'.format(np.std(psfg[0].data-psff[0].data)))
    plt.subplot(231)
    poppy.display_psf(psfg,title='Hybrid PSF',colorbar=False)
    plt.xlabel('[m]')
    plt.ylabel('[m]')
    plt.subplot(232)
    poppy.display_psf(psff,title='Fresnel PSF',colorbar=False)
    plt.xlabel('[m]')
    plt.subplot(233)
    poppy.display_psf(psfd,title='|Difference|')
    plt.xlabel('[m]')
    plt.subplot(212)
    plt.plot(gbd_pro[0],gbd_pro[1],label='Hybrid PSF Radial Average',color='red',linewidth=3,linestyle='dashdot')
    plt.plot(gbd_pro[0],fnl_pro[1],label='Fresnel PSF Radial Average',color='black',linewidth=3,linestyle='solid')
    plt.plot(gbd_pro[0],dif_pro,label='Difference',color='black',linewidth=3,linestyle='dotted')
    plt.yscale('log')
    plt.xlim([0,.0002])
    plt.ylabel('Normalized Irradiance')
    plt.xlabel('Distance [m]')
    plt.legend()
    plt.show()


FiducialCoronagraphTest()




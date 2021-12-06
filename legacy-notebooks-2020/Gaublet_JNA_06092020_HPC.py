# Package Import Section
import numpy as np
import matplotlib.pyplot as plt
import numexpr as ne
import timeit
import astropy.units as u
from scipy.special import erfc
from numba import jit

# Creating an Optical System Class To propagate rays
class GaubletOpticalSystem:

    def __init__(self,
                 epd,
                 npix,
                 dimd,
                 wavelength,
                 numbeamlets):

        basesys = np.array([[1.0,0,0,0],
                            [0,1.0,0,0],
                            [0,0,1.0,0],
                            [0,0,0,1.0]])
        self.system = basesys
        self.epd = epd
        self.npix = npix
        self.dimd = dimd# Number of Beamlets
        self.N = numbeamlets #int(np.floor(np.pi*(W*OF/(2*wo)))*4)

        # THIS BLOCK IS HARD-CODED, CHANGE FOR FINAL VERSION ##################################################

        # Beamlet Parameters
        self.wl = wavelength * u.meter # beamlet wavelength
        OF = 1.7 # Overlap Factor
        wo = 10.0*self.wl # beamlet waist
        zr = np.pi*wo**2.0/self.wl

        # Create system with a circular aperture for testing - arbitrary for now
        self.epd = epd
        self.npix = npix
        dimd = dimd

        # Create List of Positions (X,Y) in a Fibbonacci Sampled Spiral Circular Aperture
        c = np.array([0,0]) # XY offset from a spiral
        R = self.epd*np.sqrt(np.linspace(1/2,self.N-1/2,self.N))/np.sqrt(self.N-1/2)
        T = 4/(1+np.sqrt(5))*np.pi*np.linspace(1,self.N,self.N);
        X = c[0] +R*np.cos(T)
        Y = c[1] +R*np.sin(T)

        # THIS BLOCK IS HARD-CODED, CHANGE FOR FINAL VERSION ##################################################


        # Define a Q Matrix - diagonal zero for nonastigmatic, nonrotated case
        qxx = 1.0/(1j*zr)
        qxy = 0.0
        qyx = 0.0
        qyy = 1.0/(1j*zr)
        self.Q = np.array([[qxx,qxy],
                            [qyx,qyy]],dtype='complex') # Defines the matrix of inverse q parameters

        # Create the Base Rays
        self.baserays = np.array([X,
                                  Y,
                                  0*X,
                                  0*Y]) # slopes are all 0 for the base ray

    def add_optic(self,efl):
        efl = efl

        # Focusing matrix
        optic = np.array([[1.0,0,0,0],
                  [0,1.0,0,0],
                  [-1.0/float(efl),0,1.0,0],
                  [0,-1.0/float(efl),0,1.0]])
        self.system = np.matmul(optic,self.system)

    def add_distance(self,distance,index):
        distance = distance
        index = index

        # Propagation matrix
        propg = np.array([[1.0,0,float(distance)/float(index),0],
                          [0,1.0,0,float(distance)/float(index)],
                          [0,0,1.0,0],
                          [0,0,0,1.0]])
        self.system = np.matmul(propg,self.system)

    def propagate(self):

        # Propagate the base rays
        prop = np.matmul(self.system,self.baserays)

        # Propagate the Q matrix
        A = self.system[0:2,0:2]
        B = self.system[0:2,2:4]
        C = self.system[2:4,0:2]
        D = self.system[2:4,2:4]
        Qprop_n = (C + np.matmul(D,self.Q))
        Qprop_d = np.linalg.inv(A+np.matmul(B,self.Q))
        Qprop   = np.matmul(Qprop_n,Qprop_d)
        return Qprop,prop

    #def display():


class GaubletWavefront:

    def __init__(self,
                 wavelength,
                 numbeamlets,
                 npix,
                 dimension,
                 proprays,
                 baserays,
                 Qorig,
                 Qprop,
                 system):

        self.wavelength = wavelength
        self.numbeamlets = numbeamlets
        self.npix = npix
        self.dimension = dimension
        self.proprays = proprays
        self.baserays = baserays
        self.Q = Qorig
        self.Qprop = Qprop
        self.system = system
        u = np.linspace(-self.dimension,self.dimension,self.npix)
        v = np.linspace(-self.dimension,self.dimension,self.npix)
        self.u,self.v = np.meshgrid(u,v)

        # pre-define a datacube to dump the Gaublet phase in
        self.Dphase = np.zeros([npix,npix,self.numbeamlets],dtype='complex')

    def Phasecalc(self): # returns datacube of gaublet phases


        # BIG CHANGE HAPPENED HERE
        # Eikonal along propagation axis is JUST THE Z AXIS
        lo = 1*self.system[0,2] + 0*(np.sqrt((self.proprays[0,:] - self.baserays[0,:])**2 + (self.proprays[1,:] - self.baserays[1,:])**2 +self.system[0,2]**2))
        A  = self.system[0:2,0:2]
        B  = self.system[0:2,2:4]
        phase = self.Phasecube(self.numbeamlets,
                               self.dimension,
                               self.npix,
                               self.proprays,
                               self.wavelength,
                               self.Qprop,
                               lo,
                               self.Dphase,
                               self.u,
                               self.v)
        phasor = ne.evaluate('exp(phase)')
        Ephase = np.sum(phasor,axis=2)/np.sqrt(np.linalg.norm(A+np.matmul(B,self.Q)))

        return Ephase

    @staticmethod
    @jit(nopython=True)
    def Phasecube(numbeamlets,dimension,npix,proprays,wavelength,Qprop,lo,Dphase,u,v):
        for ind in range(numbeamlets):

            up = u - proprays[0,ind]
            vp = v - proprays[1,ind]

            tran_phase = ((-1j*(np.pi/wavelength))*(Qprop[0,0]*np.square(up) + (Qprop[1,0] + Qprop[0,1])*up*vp + Qprop[1,1]*np.square(vp)))
            long_phase = (-1j*(2.0*np.pi/wavelength)*lo[ind])
            Dphase[:,:,ind] = tran_phase+long_phase

        return Dphase

    def display(self,field):

        self.field = field

        # displays field amplitude, phase, and irradiance
        u = np.linspace(-self.dimension,self.dimension,self.npix)
        v = np.linspace(-self.dimension,self.dimension,self.npix)

        plt.figure(figsize=[17,9])
        plt.subplot(121)
        plt.set_cmap('gray')
        plt.pcolor(u,v,np.log(np.abs(self.field)))
        plt.title('Field Amplitude')
        plt.xlabel('Detector Dimension [m]')
        plt.ylabel('Detector Dimension [m]')

        plt.subplot(122)
        plt.pcolor(u,v,np.angle((self.field)))
        plt.title('Field Phase')
        plt.xlabel('Detector Dimension [m]')
        plt.ylabel('Detector Dimension [m]')
        plt.show()



# Test System - EPD MIGHT BE EPR
osys = GaubletOpticalSystem(epd=5e-4,npix=512,dimd=5e-4,wavelength=2.2e-6,numbeamlets=1000)
osys.add_optic(efl=50e-4)
osys.add_distance(distance=50e-4,index=1)
Qp,prop = osys.propagate()
gwfr = GaubletWavefront(wavelength=osys.wl,numbeamlets=osys.N,npix=osys.npix,dimension=osys.dimd,proprays=prop,baserays=osys.baserays,Qorig=osys.Q,Qprop=Qp,system = osys.system)
Dfield = gwfr.Phasecalc()
gwfr.display(field=Dfield)

# (2.2e-6 *(3.0/2.0))
# np.conj(self.field)
# settings produce a cool phase pattern
# epd = 5e-3
# npix = 512
# dimd = 1e-2
# wl = 2.2e-6
# numbeamlets = 3000

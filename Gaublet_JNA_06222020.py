# Package Import Section
import numpy as np
import matplotlib.pyplot as plt
import numexpr as ne
import timeit
import astropy.units as u
from scipy.special import erfc
from numba import jit
from numpy.fft import fft

samplescheme = 'fi'

# Creating an Optical System Class To propagate rays
class GaubletOpticalSystem:

    def __init__(self,
                 epd,
                 npix,
                 dimd,
                 wavelength):

        basesys = np.array([[1.0,0,0,0],
                            [0,1.0,0,0],
                            [0,0,1.0,0],
                            [0,0,0,1.0]])
        self.system = basesys
        self.epd = epd
        self.npix = npix
        self.dimd = dimd

        # THIS BLOCK IS HARD-CODED, CHANGE FOR FINAL VERSION ##################################################

        # Beamlet Parameters
        self.wl = wavelength# beamlet wavelength
        OF = 2 # Overlap Factor
        wo = 3.0*self.wl # beamlet waist
        zr = np.pi*wo**2.0/self.wl

        

        if samplescheme == 'fib':

          # Create List of Positions (X,Y) in a Fibbonacci Sampled Spiral Circular Aperture
          self.N = np.int(np.round(np.pi*((self.epd/2.0)*OF/(wo))*9.0)) # EXPERIMENTAL
          print('numbeamlets = ',self.N)
          c = np.array([0,0]) # XY offset from a spiral
          R = (self.epd/2)*np.sqrt(np.linspace(1/2,self.N-1/2,self.N))/np.sqrt(self.N-1/2) # linear space should start at 1/2
          T = 4/(1+np.sqrt(5))*np.pi*np.linspace(1,self.N,self.N);
          X = c[0] +R*np.cos(T)
          Y = c[1] +R*np.sin(T)

        else:
          # default to grid samplings
          self.N = int(round(self.epd*OF/(2*wo)))
          print('numbeamlets across grid = ',self.N)
          x = np.linspace(-self.epd/2,self.epd/2,self.N)
          y = np.linspace(-self.epd/2,self.epd/2,self.N)
          x,y = np.meshgrid(x,y)
          X = np.concatenate(x).flatten('F')
          Y = np.concatenate(y).flatten('F')
          self.N = self.N**2
          print('total numbeamlets = ',self.N)
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
                                  0.0*X,
                                  0.0*Y]) # slopes are all 0 for the base ray

    def add_optic(self,efl):
        efl = efl

        # Focusing matrix
        optic = np.array([[1.0,0.0,0.0,0.0],
                  [0.0,1.0,0.0,0.0],
                  [-1.0/float(efl),0.0,1.0,0.0],
                  [0.0,-1.0/float(efl),0.0,1.0]])
        self.system = np.matmul(optic,self.system)

    def add_distance(self,distance,index):
        distance = distance
        index = index

        # Propagation matrix
        propg = np.array([[1.0,0.0,float(distance)/float(index),0.0],
                          [0.0,1.0,0.0,float(distance)/float(index)],
                          [0.0,0.0,1.0,0.0],
                          [0.0,0.0,0.0,1.0]])
        self.system = np.matmul(propg,self.system)

    def add_aperture(self,shape,radi):

      if shape == 'square':
        print('hi')
      elif shape == 'circle':
        print('yo')

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
        #print(self.baserays[0:2])
        #print(np.linalg.inv(np.matmul(C,self.Q) + D))
        #self.P_pram = np.matmul(np.linalg.inv(np.matmul(C,self.Q) + D),self.baserays[0:2])
        self.P_pram = np.matmul(np.linalg.inv(np.matmul(C,np.linalg.inv(self.Q))+D),self.baserays[0:2])
        #np.matmul(np.linalg.inv(C+np.matmul(D,self.Q)),np.matmul(self.Q,self.baserays[0:2])) # minus prop
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
                 system,
                 P_pram):

        self.wavelength = wavelength
        self.numbeamlets = numbeamlets
        self.npix = npix
        self.dimension = dimension
        self.proprays = proprays
        self.baserays = baserays
        self.Q = Qorig
        self.Qprop = Qprop
        self.system = system
        u = np.linspace(-self.dimension/2,self.dimension/2,self.npix)
        v = np.linspace(-self.dimension/2,self.dimension/2,self.npix)
        self.u,self.v = np.meshgrid(u,v)
        self.P_pram = P_pram

        # pre-define a datacube to dump the Gaublet phase in
        self.Dphase = np.zeros([npix,npix,self.numbeamlets],dtype='complex')


    def Phasecalc(self): # returns datacube of gaublet phases

        # BIG CHANGE HAPPENED HERE
        # Eikonal along propagation axis is JUST THE Z AXIS
        lo = self.system[0,2]
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
                               self.v,
                               self.P_pram[0,:],
                               self.P_pram[1,:],
                               self.baserays)
        phasor = ne.evaluate('exp(phase)')
        Ephase = np.sum(phasor,axis=2)/np.sqrt(np.linalg.norm(A+np.matmul(B,self.Q)))

        return Ephase

    @staticmethod
    @jit(nopython=True)
    def Phasecube(numbeamlets,dimension,npix,proprays,wavelength,Qprop,lo,Dphase,u,v,P_x,P_y,baserays):
        for ind in range(numbeamlets):
            #print(P_x)
            #print(P_y)

            up = u - P_x[ind] # -  # 
            vp = v - P_y[ind] #  #
            #print(up)
            #print(vp)
            guoy_phase = -1j*np.arctan(lo/np.real(Qprop[0,0]))
            tran_phase = (-1j*(np.pi/wavelength))*(Qprop[0,0]*up**2 + (Qprop[1,0] + Qprop[0,1])*up*vp + Qprop[1,1]*vp**2)
            long_phase = -1j*(2.0*np.pi/wavelength)*lo
            Dphase[:,:,ind] = tran_phase+long_phase+guoy_phase

        return Dphase

    def display(self,field):

        x = np.linspace(-self.dimension/2,self.dimension/2,self.npix)

        self.field = field

        # displays field amplitude, phase, and irradiance
        #u = np.linspace(-self.dimension/2,self.dimension/2,self.npix)
        #v = np.linspace(-self.dimension/2,self.dimension/2,self.npix)
        #u,v = np.meshgrid(u,v)

        #print(self.u)

        plt.figure(1,figsize=[17,9])
        plt.subplot(1,2,1)
        plt.set_cmap('gray')
        plt.imshow((np.abs(self.field)))
        plt.title('Field Amplitude')
        plt.xlabel('Detector Dimension [m]')
        plt.ylabel('Detector Dimension [m]')

        plt.subplot(1,2,2)
        plt.imshow(np.angle((self.field)))
        plt.title('Field Phase')
        plt.xlabel('Detector Dimension [m]')
        plt.ylabel('Detector Dimension [m]')
        plt.show()

        plt.figure(2,figsize=[17,9])
        plt.subplot(121)
        plt.plot(x,np.abs(self.field[int(self.npix/2),:]))
        plt.title('Amplitude Cross Section X')
        plt.xlabel('Detector Dimension [m]')
        plt.ylabel('Amplitude')

        plt.subplot(1,2,2)
        plt.plot(x,np.abs((self.field[:,int(self.npix/2)])))
        plt.title('Amplitude Cross Section Y')
        plt.xlabel('Detector Dimension [m]')
        plt.ylabel('Amplitude')
        plt.show()


# Test System - EPD MIGHT BE EPR
osys = GaubletOpticalSystem(epd=2.5e-4,npix=512,dimd=5e-4,wavelength=2.2e-6) # 1e-5
osys.add_optic(efl=2.5e-4)
osys.add_distance(distance=2.2e-6,index=1)
Qp,prop = osys.propagate()
gwfr = GaubletWavefront(wavelength=osys.wl,numbeamlets=osys.N,npix=osys.npix,dimension=osys.dimd,proprays=prop,baserays=osys.baserays,Qorig=osys.Q,Qprop=Qp,system = osys.system,P_pram=osys.P_pram)
Dfield = gwfr.Phasecalc()
gwfr.display(field=Dfield)
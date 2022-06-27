# Utility
import numpy as np
import matplotlib.pyplot as plt

# Packages to Test
import raytrace as rt
import gbd as gb
import differential as dif

# raytrace tests

# gbd tests

# The test case should work for a non-orthogonal optical system
npix = 256
wl = 1e-6 # meters
wo = 100000*wl # meters
d = 100
efl = 100
zr = np.pi * wo**2 / wl # meters
q  = 1j*zr
qinv = q**-1
k = 2*np.pi/wl

# an orthogonal beamlet
Qinv = np.array([[qinv,0],
                 [0,qinv]])

# non-orthogonal optical system
ABCD = np.array([[1e-20,0,d,0],
                 [0,1e-20,0,d],
                 [-1/efl,0,1,0],
                 [0,-1/efl,0,1]])

A = ABCD[0:2,0:2]
B = ABCD[0:2,2:4]
C = ABCD[2:4,0:2]
D = ABCD[2:4,2:4]

# decenter parameter
rd = np.array([1/2,1/2])
xd = rd[0]
yd = rd[1]
ray_in = np.array([xd,yd,0,0])
ray_out = ABCD @ ray_in

# Original Beamlet
x1 = np.linspace(-1,1,npix)
x1,y1 = np.meshgrid(x1,x1)

x2 = np.linspace(-1e-2,1e-2,npix)
x2,y2 = np.meshgrid(x2,x2)

x2 -= ray_out[0]
y2 -= ray_out[1]

x1 -= xd
y1 -= yd

Qpinv = (C + D @ Qinv) @ np.linalg.inv(A + B @ Qinv)
demat = np.linalg.inv(np.linalg.inv(Qinv) + np.linalg.inv(A) @ B)
crmat = np.linalg.inv(A @ np.linalg.inv(Qinv) + B)


E1 = np.exp(-1j*k/2 * (x1**2 * Qinv[0,0] + (Qinv[1,0] + Qinv[0,1])*x1*y1 + Qinv[1,1]*y1**2))

E2  = np.exp(-1j*k/2 * (x2**2 * Qpinv[0,0] + (Qpinv[1,0] + Qpinv[0,1])*x2*y2 + Qpinv[1,1]*y2**2))
E2 *= np.exp(-1j*k/2 * (xd**2 * demat[0,0] + (demat[1,0] + demat[0,1])*xd*yd + demat[1,1]*yd**2))
E2 *= np.exp(1j*k * (xd*crmat[0,0]*x2 + xd*crmat[0,1]*y2 + yd*crmat[1,0]*x2 + yd*crmat[1,1]*y2))
E2 *= np.sqrt(np.linalg.det(A + B @ Qinv))

plt.figure()
plt.imshow(np.abs(E1*np.conj(E1)))
plt.title('E1 Irradiance')
plt.colorbar()
plt.show()

plt.figure()
plt.imshow(np.abs(E2*np.conj(E2)))
plt.title('E2 Irradiance')
plt.colorbar()
plt.show()












#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 14:51:18 2022

@author: jashcraft
"""

import numpy as np
import numpy
import matplotlib.pyplot as plt
import gfuncs as gf

npix = 512
scale = 4e-1
size = 25.4e-3
erx = []
ery = []

# do percent error calculation for the first 37 Zernikes

# x = (np.linspace(-size/2,size/2,npix) + 1e-20)
# x,y = np.meshgrid(x,x)
def zernike_j(index,size):
    x = np.linspace(-1,1,size) + 1e-20
    x,y = np.meshgrid(x,x)
    
    if index == 0:
        z = 0*x+1
    
    elif index == 1:
        z = 2*y
    
    elif index == 2:
        z = 2*x
        
    elif index == 3:
        z = np.sqrt(6)*2*y*x
        
    elif index == 4:
        z = np.sqrt(3)*(2*x**2 + 2*y**2 - 1)
        
    elif index == 5:
        z = np.sqrt(6)*(x**2 - y**2)
        
    elif index == 6:
        z = np.sqrt(8)*(3*x**2*y - y**3)
        
    elif index == 7:
        z = np.sqrt(8)*(3*x**2*y + 3*y**3 - 2*y)
        
    elif index == 8:
        z = np.sqrt(8)*(3*x**3 + 3*x*y**2 - 2*x)
        
    elif index == 9:
        z = np.sqrt(8)*(x**3 - 3*x*y**2)
        
    elif index == 10:
        z = np.sqrt(10)*(4*x**3*y - 4*x*y**3)
        
    elif index == 11:
        z = np.sqrt(10)*(8*x**3*y + 8*x*y**3 - 6*x*y)
        
    elif index == 12:
        z = np.sqrt(5)*(6*x**4 + 12*x**2*y**2 + 6*y**4 - 6*x**2 - 6*y**2 + 1)
        
    return z
    
    

def zernike_derivative_j(index,size):
    x = np.linspace(-1,1,size) + 1e-20
    x,y = np.meshgrid(x,x)
    
    if index == 0:
        dx = 0*y
        dy = 0*x
    
    elif index == 1:
        dx = 0*y
        dy = 2+0*x
    
    elif index == 2:
        dx = 2+0*y
        dy = 0*x
        
    elif index == 3:
        dx = np.sqrt(6)*2*y
        dy = np.sqrt(6)*2*x
        
    elif index == 4:
        dx = np.sqrt(3)*4*x
        dy = np.sqrt(3)*4*y
        
    elif index == 5:
        dx = np.sqrt(6)*2*x
        dy = -np.sqrt(6)*2*y
        
    elif index == 6:
        dx = np.sqrt(8)*6*x*y
        dy = np.sqrt(8)*(3*x**2 - 3*y**2)
        
    elif index == 7:
        dx = np.sqrt(8)*6*x*y
        dy = np.sqrt(8)*(3*x**2 + 9*y**2 - 2)
        
    elif index == 8:
        dx = np.sqrt(8)*(9*x**2 + 3*y**2 - 2)
        dy = np.sqrt(8)*6*x*y
        
    elif index == 9:
        dx = np.sqrt(8)*(3*x**2 - 3*y**2)
        dy = -np.sqrt(8)*6*x*y
        
    elif index == 10:
        dx = np.sqrt(10)*(12*x**2*y - 4*y**3)
        dy = np.sqrt(10)*(4*x**3 - 12*x*y**2)
        
    elif index == 11:
        dx = np.sqrt(10)*(24*x**2*y + 8*y**3 - 6*y)
        dy = np.sqrt(10)*(8*x**3 + 24*x*y**2 - 6*x)
        
    elif index == 12:
        dx = np.sqrt(5)*(24*x**3 + 24*x*y**2 - 12*x)
        dy = np.sqrt(5)*(24*x**2*y + 24*y**3 - 12*y)
        
    return dx,dy

# for i in range(1,11):
    
#     # compute the zernike derivative analytically
#     dx_a,dy_a = zernike_derivative_j(i,npix)
    
#     # compute the zernike derivative computationally
#     zern = gf.zernike_j(i, npix)
#     dy_c,dx_c = np.gradient(zern,2/npix)
    
#     # compute percent error
#     erx.append(np.mean(100*(dx_c-dx_a)/dx_a))
#     ery.append(np.mean(100*(dy_c-dy_a)/dy_a))


#z5  = np.sqrt(6)*((x)**2 - (y)**2) 
#dxz5_analytical = np.sqrt(6)*2*(x)
#dyz5_analytical = -np.sqrt(6)*2*(y)
#dyz5_computed,dxz5_computed = np.gradient(z5,size/(npix))

# different shapes at Z8

index = 7
vmin = 1e-1
vmax = 3e-1

z = zernike_j(index, npix)
dy_c,dx_c = np.gradient(z,2/npix)
dx_a,dy_a = zernike_derivative_j(index,npix)

perror_x = 100*(dx_c-dx_a)/dx_a 
perror_y = 100*(dy_c-dy_a)/dy_a 

dy_c += 1e-20
dx_c += 1e-20
dy_a += 1e-20
dx_a += 1e-20

plt.figure()
plt.title('Zernike {} Under Test'.format(index))
plt.imshow(z)
plt.colorbar()
plt.show()

plt.figure(figsize=[16,7])
plt.subplot(121)
plt.imshow(dx_a,vmin=-scale*2,vmax=scale*2)
plt.title('Derivative in X - derived')
plt.colorbar()
plt.subplot(122)
plt.imshow(dy_a,vmin=-scale*2,vmax=scale*2)
plt.title('Derivative in Y - derived')
plt.colorbar()
plt.show()

plt.figure(figsize=[16,7])
plt.subplot(121)
plt.imshow(dx_c,vmin=-scale*2,vmax=scale*2)
plt.title('Derivative in X - computed')
plt.colorbar()
plt.subplot(122)
plt.imshow(dy_c,vmin=-scale*2,vmax=scale*2)
plt.title('Derivative in Y - computed')
plt.colorbar()
plt.show()

plt.figure(figsize=[16,7])
plt.suptitle('Percent Error in np.gradient()')
plt.subplot(121)
plt.imshow(perror_x,vmin=vmin,vmax=vmax)
plt.title('dW/dx')
plt.colorbar()
plt.subplot(122)
plt.imshow(perror_y,vmin=vmin,vmax=vmax)
plt.title('dW/dy')
plt.colorbar()
plt.show()
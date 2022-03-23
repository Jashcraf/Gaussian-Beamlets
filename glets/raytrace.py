#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 11:10:36 2022

@author: jashcraft
"""
import numpy as np

def makerays(size,numrays,circle=True):
    """
    

    Parameters
    ----------
    size : float
        dimension across a square grid of rays.
    numrays : int
        number of rays across size.
    circle : bool, optional
        Whether to trace a circle of rays instead of a square. The default is True.

    Returns
    -------
    ndarray
        4xnrays array of rays.

    """
    
    # Define lists of XY coordinate pairs for square grid
    x = np.linspace(-size/2,size/2,numrays) + 1e-12
    y = np.linspace(-size/2,size/2,numrays) + 1e-12
    x,y = np.meshgrid(x,y)
    X = np.ravel(x)
    Y = np.ravel(y)
    if circle == True:
        r = np.sqrt(X**2 + Y**2)
        X = X[r <= size/2]
        Y = Y[r <= size/2]
        
    
    return np.array([X,Y,0*X,0*Y])

def IdentityMat(nrays):
    """
    
    Parameters
    ----------
    nrays : int
        number of rays across a square grid.

    Returns
    -------
    box : ndarray
        4x4xnrays array of ray transfer matrices.

    """
    
    box = np.zeros([4,4,nrays])
    
    # construct identity matrix
    box[0,0,:] = 1.
    box[1,1,:] = 1.
    box[2,2,:] = 1.
    box[3,3,:] = 1.
    
    return box
    

def AnamorphicLens(eflx,efly,nrays):
    """

    Parameters
    ----------
    eflx : float
        focal length in the x direction.
    efly : float
        focal length in the y direction.
    nrays : int
        number of rays across a square grid.

    Returns
    -------
    analens : ndarray
        4x4xnrays array of ray transfer matrices.

    """
    
    analens = IdentityMat(nrays)
    analens[2,0,:] = -1/eflx
    analens[3,1,:] = -1/efly
    
    return analens

def ThinLens(efl,nrays):
    """

    Parameters
    ----------
    efl : float
        effective focal length of the lens.
    nrays : int
        number of rays across a square grid.

    Returns
    -------
    ndarray
        4x4xnrays array of ray transfer matrices.

    """
    
    return AnamorphicLens(efl,efl,nrays)

def FreeSpace(distance,nrays):
    """

    Parameters
    ----------
    distance : float
        distance to propagate the rays.
    nrays : int
        number of rays across a square grid.

    Returns
    -------
    dmat : ndarray
        4x4xnrays array of ray transfer matrices.

    """
    
    dmat = IdentityMat(nrays)
    dmat[0,2,:] = distance
    dmat[1,3,:] = distance
    
    return dmat

def ZernikeWFE(nrays,size,rays,index,scale,plotderivs=False):
    """
    EXPERIMENTAL: Only alters the position of rays, cannot be multiplied by a 
    ray transfer matrix. Can only do the first 11 Zernike polynomials. 
    
    Derivatives of Zernike polynomials taken from "Description of Zernike Polynomials"
    
    Source: https://wp.optics.arizona.edu/visualopticslab/wp-content/uploads/sites/52/2016/08/Zernike-Notes-15Jan2016.pdf
    Author: Dr. Jim Schwiegerling

    Parameters
    ----------
    nrays : int
        number of rays across a square grid.
    size : float
        dimension across a square grid.
    rays : ndarray
        4xnrays array of rays.
    index : int
        zernike polynomial index.
    scale : float
        multiple to control the amount of zernike used.
    plotderivs : boolean, optional
        Option of whether to plot the computed derivatives. The default is False.

    Returns
    -------
    rays : ndarray
        4xnrays array of rays.

    """
    x = np.linspace(-size/2,size/2,nrays) + 1e-20 # normalize the size
    x,y = np.meshgrid(x,x)
    
    #x = x/(size/2)
    #y = y/(size/2)
    
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
        
    dx *= scale
    dy *= scale
    
    if plotderivs:
        import matplotlib.pyplot as plt
        plt.figure(figsize=[16,7])
        plt.subplot(121)
        plt.imshow(dx,vmin=-scale*10,vmax=scale*10)
        plt.title('Derivative in X - derived')
        plt.colorbar()
        plt.subplot(122)
        plt.imshow(dy,vmin=-scale*10,vmax=scale*10)
        plt.title('Derivative in Y - derived')
        plt.colorbar()
        plt.show()
    
    rays[2,:] -= np.ravel(dx)
    rays[3,:] -= np.ravel(dy)
    
    return rays

def ArbitraryWFE(nrays,size,rays,scale,zern_index=None,array=None):
    wlen = 2.2e-6
    x = np.linspace(-size/2,size/2,nrays) + 1e-20
    x,y = np.meshgrid(x,x)
    import matplotlib.pyplot as plt
    
    # This is a test case w/ know structure
    if zern_index != None:
        index = zern_index
        
        import matplotlib.pyplot as plt
        
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
            
    if array[0,0] != None:
        
        # Kind of a hacked fix for now
        z = array*scale #np.exp(1j*2*np.pi/(wlen)*array)
            
        
    #zern = zernike_j(zern_index,nrays)
    #dy,dx = np.gradient(np.angle(z)*wlen/(2*np.pi),size/(nrays),edge_order=1)
    print('Computing Derivative for Pixelscale of = ',size/nrays)
#     dy,dx = np.arctan(np.gradient(z,size/(nrays),edge_order=1))
    dy,dx = (np.gradient(z,size/(nrays),edge_order=2))


    #dx *= scale
    #dy *= scale

#     plt.figure()
#     plt.imshow(np.angle(z))
#     plt.colorbar()
#     plt.title('rms wfe applied [radians]')
#     plt.show()

#     plt.figure(figsize=[16,7])
#     plt.subplot(121)
#     plt.imshow(dx)
#     plt.title('Derivative in X')
#     plt.colorbar()
#     plt.subplot(122)
#     plt.imshow(dy)
#     plt.title('Derivative in Y')
#     plt.colorbar()
#     plt.show()

    rays[2,:] = rays[2,:]-np.ravel(dx)
    rays[3,:] = rays[3,:]-np.ravel(dy)
        
    return rays

def matmultlist(mat1,mat2):
    
    # this works for 4x4xn dimension arrays, where it multiplies the matrices element-wise by the final axis
    
    box = np.zeros(mat1.shape)
    
    for mind in range(mat1[0,0,:].size):
        
        box[:,:,mind] = np.matmul(mat1[:,:,mind],mat2[:,:,mind])
        
    return box

def matmultrays(mat1,rays):
    
    # this works for 4x4xn dimension arrays, where it multiplies the matrices element-wise by the final axis
    
    box = np.zeros(rays.shape)
    
    for mind in range(mat1[0,0,:].size):
        
        box[:,mind] = np.matmul(mat1[:,:,mind],rays[:,mind])
        
    return box

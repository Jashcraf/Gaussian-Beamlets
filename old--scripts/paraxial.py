"""
Engine for paraxial ray tracing, ray generation & ABCD matrix elements
"""
import numpy as np 

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

def matmultlist(mat1,mat2):
    
    # this works for 4x4xn dimension arrays, where it multiplies the matrices element-wise by the final axis
    
    box = np.zeros(mat1.shape)
    
    for mind in range(mat1[0,0,:].size):
        
        box[:,:,mind] = np.matmul(mat1[:,:,mind],mat2[:,:,mind])
        
    return box

def proprays(rays,system):
    
    box = np.zeros(rays.shape)
    
    for rind in range(rays[0,:].size):
        
        box[:,rind] = np.matmul(system[:,:,rind],rays[:,rind])
        
    return box
    


"""
Engine for interpretation of ray data as ray transfer matrices
"""
import numpy as np
# formatted_rays = [x,y,thx,thy],[xp,yp,thxp,thyp]

# Interpreter for given ray format, begin with OpticStudio
def readrays(data):
    
    """
    data has format [xfile,yfile,afile,bfile]
    x,y = coordinates
    a,b = direction cosines in those directions, for small angles these correspond to the z slopes
    """
    
    x = np.loadtxt(data[0])
    y = np.loadtxt(data[1])
    z = np.loadtxt(data[2])
    
    
    a = np.loadtxt(data[3])
    b = np.loadtxt(data[4])
    g = np.loadtxt(data[5])
    
    #vec_xz = np.array([a,0,g])*np.cos(g)
    #res = np.dot(vec,np.array([1,0,1]))
    
    slopex = a#np.cos(np.pi/2 - np.arccos(a))
    slopey = b#np.cos(np.pi/2 - np.arccos(b)) #np.pi/2 - np.arccos(b)
    
    
    
    # Just doing over z doesn't work, you need the distance the ray traveled - 
    # Do you need to do multiple raytraces to compute the slope? Or is it possible from direction cosines?
    
    formatted_rays = np.array([x,y,slopex,slopey])
    
    return formatted_rays

# Needs 2 raysets, before and after propagation
# In order for the data to be valued there must be a differential in both position AND angle
# Will also need a way of sorting nan values
# Also need to compute rays after the entrance pupil and after propagation, so add that to code

def orthogonal_abcd(r1,r1p,th1,th1p,r2,r2p,th2,th2p):
    
    # x case
    # Ax = ((r1p/th1) - (r2p/th2)) / ((r1/th1) - (r2/th2))
    # Ax = (r1p - (th1*r1p/th2))/(r1 - (th1*r2/th2))
    # Bx = ((r1p/r1) - (r2p/r2)) / ((th1/r1) - (th2/r2))
    # Cx = ((th1p/th1) - (th2p/th2)) / ((r1/th1) - (r2/th2))
    # Cx = (th1p - (th1*th2p/th2))/(r1 - th1*r2/th2)
    # Dx = ((th1p/r1) - (th2p/r2)) / ((th1/r1) - (th2/r2))
    
    A = (th2*r1p - th1*r2p)/(th2*r1 - th1*r2)
    B = (r1p*r2 - r1*r2p)/(th1*r2 - th2*r1)
    C = (th2*th1p - th1*th2p)/(r1*th2 - r2*th1)
    D = (r2*th1p - r1*th2p)/(r2*th1 - r1*th2)
    
    print('Check C')
    print(r1[np.isnan(C)])
    print(r2[np.isnan(C)])
    print('Check D')
    print(r1[np.isnan(D)])
    print(r2[np.isnan(D)])
    
    return A,B,C,D

def compute_orthogonal_abcd(ray1in,ray1out,ray2in,ray2out):
    
    # Needs special case for theta1 = 0 where division by zero isn't cool
    
    # grab ray coordinates
    r1p  = ray1out[0,:]
    r1   = ray1in[0,:]
    th1p = ray1out[2,:]
    th1  = ray1in[2,:]
    
    r2p  = ray2out[0,:]
    r2   = ray2in[0,:]
    th2p = ray2out[2,:]
    th2  = ray2in[2,:]
    
    Ax,Bx,Cx,Dx = orthogonal_abcd(r1,r1p,th1,th1p,r2,r2p,th2,th2p)
    
    # grab other ray coordinates
    r1p  = ray1out[1,:]
    r1   = ray1in[1,:]
    th1p = ray1out[3,:]
    th1  = ray1in[3,:]
    
    r2p  = ray2out[1,:]
    r2   = ray2in[1,:]
    th2p = ray2out[3,:]
    th2  = ray2in[3,:]
    
    Ay,By,Cy,Dy = orthogonal_abcd(r1,r1p,th1,th1p,r2,r2p,th2,th2p)
    
    return Ax,Bx,Cx,Dx,Ay,By,Cy,Dy

def compute_full_rtm(ray1in,ray1out,
                    ray2in,ray2out,
                    ray3in,ray3out,
                    ray4in,ray4out):
    
    x1 = ray1in[0,:]
    x1p = ray1out[0,:]
    u1 = ray1in[2,:]
    u1p = ray1out[2,:]
    
    y1 = ray1in[1,:]
    y1p = ray1out[1,:]
    v1 = ray1in[3,:]
    v1p = ray1out[3,:]
    
    x2 = ray2in[0,:]
    x2p = ray2out[0,:]
    u2 = ray2in[2,:]
    u2p = ray2out[2,:]
    
    y2 = ray2in[1,:]
    y2p = ray2out[1,:]
    v2 = ray2in[3,:]
    v2p = ray2out[3,:]
    
    x3 = ray3in[0,:]
    x3p = ray3out[0,:]
    u3 = ray3in[2,:]
    u3p = ray3out[2,:]
    
    y3 = ray3in[1,:]
    y3p = ray3out[1,:]
    v3 = ray3in[3,:]
    v3p = ray3out[3,:]
    
    x4 = ray4in[0,:]
    x4p = ray4out[0,:]
    u4 = ray4in[2,:]
    u4p = ray4out[2,:]
    
    y4 = ray4in[1,:]
    y4p = ray4out[1,:]
    v4 = ray4in[3,:]
    v4p = ray4out[3,:]
    
    
    Rin = np.array([[x1,x2,x3,x4],
                    [y1,y2,y3,y4],
                    [u1,u2,u3,u4],
                    [v1,v2,v3,v4]])
    
    Rout = np.array([[x1p,y1p,u1p,v1p],
                     [x2p,y2p,u2p,v2p],
                     [x3p,y3p,u3p,v3p],
                     [x4p,y4p,u4p,v4p]])
    
    if np.linalg.matrix_rank(Rin) <4:
        print('Rank of matrix too low, change differentials')
    
    Rinv = np.linalg.inv(Rin)
    Osys = np.matmul(Rout,Rinv)
    
    return Osys
    
    

def compute_full_offdiagonal_abcd(ray1in,ray1out,
                                  ray2in,ray2out,
                                  ray3in,ray3out,
                                  ray4in,ray4out):
    
    x1 = ray1in[0,:]
    x1p = ray1out[0,:]
    u1 = ray1in[2,:]
    u1p = ray1out[2,:]
    
    y1 = ray1in[1,:]
    y1p = ray1out[1,:]
    v1 = ray1in[3,:]
    v1p = ray1out[3,:]
    
    x2 = ray2in[0,:]
    x2p = ray2out[0,:]
    u2 = ray2in[2,:]
    u2p = ray2out[2,:]
    
    y2 = ray2in[1,:]
    y2p = ray2out[1,:]
    v2 = ray2in[3,:]
    v2p = ray2out[3,:]
    
    x3 = ray3in[0,:]
    x3p = ray3out[0,:]
    u3 = ray3in[2,:]
    u3p = ray3out[2,:]
    
    y3 = ray3in[1,:]
    y3p = ray3out[1,:]
    v3 = ray3in[3,:]
    v3p = ray3out[3,:]
    
    x4 = ray4in[0,:]
    x4p = ray4out[0,:]
    u4 = ray4in[2,:]
    u4p = ray4out[2,:]
    
    y4 = ray4in[1,:]
    y4p = ray4out[1,:]
    v4 = ray4in[3,:]
    v4p = ray4out[3,:]
    
    # Python Solutions Sympy
    Axx = (u1*v2*x3p*y4 - u1*v2*x4p*y3 - u1*v3*x2p*y4 + u1*v3*x4p*y2 + u1*v4*x2p*y3 - u1*v4*x3p*y2 - u2*v1*x3p*y4 + u2*v1*x4p*y3 + u2*v3*x1p*y4 - u2*v3*x4p*y1 - u2*v4*x1p*y3 + u2*v4*x3p*y1 + u3*v1*x2p*y4 - u3*v1*x4p*y2 - u3*v2*x1p*y4 + u3*v2*x4p*y1 + u3*v4*x1p*y2 - u3*v4*x2p*y1 - u4*v1*x2p*y3 + u4*v1*x3p*y2 + u4*v2*x1p*y3 - u4*v2*x3p*y1 - u4*v3*x1p*y2 + u4*v3*x2p*y1)/(u1*v2*x3*y4 - u1*v2*x4*y3 - u1*v3*x2*y4 + u1*v3*x4*y2 + u1*v4*x2*y3 - u1*v4*x3*y2 - u2*v1*x3*y4 + u2*v1*x4*y3 + u2*v3*x1*y4 - u2*v3*x4*y1 - u2*v4*x1*y3 + u2*v4*x3*y1 + u3*v1*x2*y4 - u3*v1*x4*y2 - u3*v2*x1*y4 + u3*v2*x4*y1 + u3*v4*x1*y2 - u3*v4*x2*y1 - u4*v1*x2*y3 + u4*v1*x3*y2 + u4*v2*x1*y3 - u4*v2*x3*y1 - u4*v3*x1*y2 + u4*v3*x2*y1)
    Axy = (u1*v2*x3*x4p - u1*v2*x3p*x4 - u1*v3*x2*x4p + u1*v3*x2p*x4 + u1*v4*x2*x3p - u1*v4*x2p*x3 - u2*v1*x3*x4p + u2*v1*x3p*x4 + u2*v3*x1*x4p - u2*v3*x1p*x4 - u2*v4*x1*x3p + u2*v4*x1p*x3 + u3*v1*x2*x4p - u3*v1*x2p*x4 - u3*v2*x1*x4p + u3*v2*x1p*x4 + u3*v4*x1*x2p - u3*v4*x1p*x2 - u4*v1*x2*x3p + u4*v1*x2p*x3 + u4*v2*x1*x3p - u4*v2*x1p*x3 - u4*v3*x1*x2p + u4*v3*x1p*x2)/(u1*v2*x3*y4 - u1*v2*x4*y3 - u1*v3*x2*y4 + u1*v3*x4*y2 + u1*v4*x2*y3 - u1*v4*x3*y2 - u2*v1*x3*y4 + u2*v1*x4*y3 + u2*v3*x1*y4 - u2*v3*x4*y1 - u2*v4*x1*y3 + u2*v4*x3*y1 + u3*v1*x2*y4 - u3*v1*x4*y2 - u3*v2*x1*y4 + u3*v2*x4*y1 + u3*v4*x1*y2 - u3*v4*x2*y1 - u4*v1*x2*y3 + u4*v1*x3*y2 + u4*v2*x1*y3 - u4*v2*x3*y1 - u4*v3*x1*y2 + u4*v3*x2*y1)
    Bxx = (v1*x2*x3p*y4 - v1*x2*x4p*y3 - v1*x2p*x3*y4 + v1*x2p*x4*y3 + v1*x3*x4p*y2 - v1*x3p*x4*y2 - v2*x1*x3p*y4 + v2*x1*x4p*y3 + v2*x1p*x3*y4 - v2*x1p*x4*y3 - v2*x3*x4p*y1 + v2*x3p*x4*y1 + v3*x1*x2p*y4 - v3*x1*x4p*y2 - v3*x1p*x2*y4 + v3*x1p*x4*y2 + v3*x2*x4p*y1 - v3*x2p*x4*y1 - v4*x1*x2p*y3 + v4*x1*x3p*y2 + v4*x1p*x2*y3 - v4*x1p*x3*y2 - v4*x2*x3p*y1 + v4*x2p*x3*y1)/(u1*v2*x3*y4 - u1*v2*x4*y3 - u1*v3*x2*y4 + u1*v3*x4*y2 + u1*v4*x2*y3 - u1*v4*x3*y2 - u2*v1*x3*y4 + u2*v1*x4*y3 + u2*v3*x1*y4 - u2*v3*x4*y1 - u2*v4*x1*y3 + u2*v4*x3*y1 + u3*v1*x2*y4 - u3*v1*x4*y2 - u3*v2*x1*y4 + u3*v2*x4*y1 + u3*v4*x1*y2 - u3*v4*x2*y1 - u4*v1*x2*y3 + u4*v1*x3*y2 + u4*v2*x1*y3 - u4*v2*x3*y1 - u4*v3*x1*y2 + u4*v3*x2*y1)
    Bxy = (-u1*x2*x3p*y4 + u1*x2*x4p*y3 + u1*x2p*x3*y4 - u1*x2p*x4*y3 - u1*x3*x4p*y2 + u1*x3p*x4*y2 + u2*x1*x3p*y4 - u2*x1*x4p*y3 - u2*x1p*x3*y4 + u2*x1p*x4*y3 + u2*x3*x4p*y1 - u2*x3p*x4*y1 - u3*x1*x2p*y4 + u3*x1*x4p*y2 + u3*x1p*x2*y4 - u3*x1p*x4*y2 - u3*x2*x4p*y1 + u3*x2p*x4*y1 + u4*x1*x2p*y3 - u4*x1*x3p*y2 - u4*x1p*x2*y3 + u4*x1p*x3*y2 + u4*x2*x3p*y1 - u4*x2p*x3*y1)/(u1*v2*x3*y4 - u1*v2*x4*y3 - u1*v3*x2*y4 + u1*v3*x4*y2 + u1*v4*x2*y3 - u1*v4*x3*y2 - u2*v1*x3*y4 + u2*v1*x4*y3 + u2*v3*x1*y4 - u2*v3*x4*y1 - u2*v4*x1*y3 + u2*v4*x3*y1 + u3*v1*x2*y4 - u3*v1*x4*y2 - u3*v2*x1*y4 + u3*v2*x4*y1 + u3*v4*x1*y2 - u3*v4*x2*y1 - u4*v1*x2*y3 + u4*v1*x3*y2 + u4*v2*x1*y3 - u4*v2*x3*y1 - u4*v3*x1*y2 + u4*v3*x2*y1)
    
    Ayx = (-u1*v2*y3*y4p + u1*v2*y3p*y4 + u1*v3*y2*y4p - u1*v3*y2p*y4 - u1*v4*y2*y3p + u1*v4*y2p*y3 + u2*v1*y3*y4p - u2*v1*y3p*y4 - u2*v3*y1*y4p + u2*v3*y1p*y4 + u2*v4*y1*y3p - u2*v4*y1p*y3 - u3*v1*y2*y4p + u3*v1*y2p*y4 + u3*v2*y1*y4p - u3*v2*y1p*y4 - u3*v4*y1*y2p + u3*v4*y1p*y2 + u4*v1*y2*y3p - u4*v1*y2p*y3 - u4*v2*y1*y3p + u4*v2*y1p*y3 + u4*v3*y1*y2p - u4*v3*y1p*y2)/(u1*v2*x3*y4 - u1*v2*x4*y3 - u1*v3*x2*y4 + u1*v3*x4*y2 + u1*v4*x2*y3 - u1*v4*x3*y2 - u2*v1*x3*y4 + u2*v1*x4*y3 + u2*v3*x1*y4 - u2*v3*x4*y1 - u2*v4*x1*y3 + u2*v4*x3*y1 + u3*v1*x2*y4 - u3*v1*x4*y2 - u3*v2*x1*y4 + u3*v2*x4*y1 + u3*v4*x1*y2 - u3*v4*x2*y1 - u4*v1*x2*y3 + u4*v1*x3*y2 + u4*v2*x1*y3 - u4*v2*x3*y1 - u4*v3*x1*y2 + u4*v3*x2*y1)
    Ayy = (u1*v2*x3*y4p - u1*v2*x4*y3p - u1*v3*x2*y4p + u1*v3*x4*y2p + u1*v4*x2*y3p - u1*v4*x3*y2p - u2*v1*x3*y4p + u2*v1*x4*y3p + u2*v3*x1*y4p - u2*v3*x4*y1p - u2*v4*x1*y3p + u2*v4*x3*y1p + u3*v1*x2*y4p - u3*v1*x4*y2p - u3*v2*x1*y4p + u3*v2*x4*y1p + u3*v4*x1*y2p - u3*v4*x2*y1p - u4*v1*x2*y3p + u4*v1*x3*y2p + u4*v2*x1*y3p - u4*v2*x3*y1p - u4*v3*x1*y2p + u4*v3*x2*y1p)/(u1*v2*x3*y4 - u1*v2*x4*y3 - u1*v3*x2*y4 + u1*v3*x4*y2 + u1*v4*x2*y3 - u1*v4*x3*y2 - u2*v1*x3*y4 + u2*v1*x4*y3 + u2*v3*x1*y4 - u2*v3*x4*y1 - u2*v4*x1*y3 + u2*v4*x3*y1 + u3*v1*x2*y4 - u3*v1*x4*y2 - u3*v2*x1*y4 + u3*v2*x4*y1 + u3*v4*x1*y2 - u3*v4*x2*y1 - u4*v1*x2*y3 + u4*v1*x3*y2 + u4*v2*x1*y3 - u4*v2*x3*y1 - u4*v3*x1*y2 + u4*v3*x2*y1)
    Byx = (-v1*x2*y3*y4p + v1*x2*y3p*y4 + v1*x3*y2*y4p - v1*x3*y2p*y4 - v1*x4*y2*y3p + v1*x4*y2p*y3 + v2*x1*y3*y4p - v2*x1*y3p*y4 - v2*x3*y1*y4p + v2*x3*y1p*y4 + v2*x4*y1*y3p - v2*x4*y1p*y3 - v3*x1*y2*y4p + v3*x1*y2p*y4 + v3*x2*y1*y4p - v3*x2*y1p*y4 - v3*x4*y1*y2p + v3*x4*y1p*y2 + v4*x1*y2*y3p - v4*x1*y2p*y3 - v4*x2*y1*y3p + v4*x2*y1p*y3 + v4*x3*y1*y2p - v4*x3*y1p*y2)/(u1*v2*x3*y4 - u1*v2*x4*y3 - u1*v3*x2*y4 + u1*v3*x4*y2 + u1*v4*x2*y3 - u1*v4*x3*y2 - u2*v1*x3*y4 + u2*v1*x4*y3 + u2*v3*x1*y4 - u2*v3*x4*y1 - u2*v4*x1*y3 + u2*v4*x3*y1 + u3*v1*x2*y4 - u3*v1*x4*y2 - u3*v2*x1*y4 + u3*v2*x4*y1 + u3*v4*x1*y2 - u3*v4*x2*y1 - u4*v1*x2*y3 + u4*v1*x3*y2 + u4*v2*x1*y3 - u4*v2*x3*y1 - u4*v3*x1*y2 + u4*v3*x2*y1)
    Byy = (u1*x2*y3*y4p - u1*x2*y3p*y4 - u1*x3*y2*y4p + u1*x3*y2p*y4 + u1*x4*y2*y3p - u1*x4*y2p*y3 - u2*x1*y3*y4p + u2*x1*y3p*y4 + u2*x3*y1*y4p - u2*x3*y1p*y4 - u2*x4*y1*y3p + u2*x4*y1p*y3 + u3*x1*y2*y4p - u3*x1*y2p*y4 - u3*x2*y1*y4p + u3*x2*y1p*y4 + u3*x4*y1*y2p - u3*x4*y1p*y2 - u4*x1*y2*y3p + u4*x1*y2p*y3 + u4*x2*y1*y3p - u4*x2*y1p*y3 - u4*x3*y1*y2p + u4*x3*y1p*y2)/(u1*v2*x3*y4 - u1*v2*x4*y3 - u1*v3*x2*y4 + u1*v3*x4*y2 + u1*v4*x2*y3 - u1*v4*x3*y2 - u2*v1*x3*y4 + u2*v1*x4*y3 + u2*v3*x1*y4 - u2*v3*x4*y1 - u2*v4*x1*y3 + u2*v4*x3*y1 + u3*v1*x2*y4 - u3*v1*x4*y2 - u3*v2*x1*y4 + u3*v2*x4*y1 + u3*v4*x1*y2 - u3*v4*x2*y1 - u4*v1*x2*y3 + u4*v1*x3*y2 + u4*v2*x1*y3 - u4*v2*x3*y1 - u4*v3*x1*y2 + u4*v3*x2*y1)
    
    Cxx = (-u1*u2p*v3*y4 + u1*u2p*v4*y3 + u1*u3p*v2*y4 - u1*u3p*v4*y2 - u1*u4p*v2*y3 + u1*u4p*v3*y2 + u1p*u2*v3*y4 - u1p*u2*v4*y3 - u1p*u3*v2*y4 + u1p*u3*v4*y2 + u1p*u4*v2*y3 - u1p*u4*v3*y2 - u2*u3p*v1*y4 + u2*u3p*v4*y1 + u2*u4p*v1*y3 - u2*u4p*v3*y1 + u2p*u3*v1*y4 - u2p*u3*v4*y1 - u2p*u4*v1*y3 + u2p*u4*v3*y1 - u3*u4p*v1*y2 + u3*u4p*v2*y1 + u3p*u4*v1*y2 - u3p*u4*v2*y1)/(u1*v2*x3*y4 - u1*v2*x4*y3 - u1*v3*x2*y4 + u1*v3*x4*y2 + u1*v4*x2*y3 - u1*v4*x3*y2 - u2*v1*x3*y4 + u2*v1*x4*y3 + u2*v3*x1*y4 - u2*v3*x4*y1 - u2*v4*x1*y3 + u2*v4*x3*y1 + u3*v1*x2*y4 - u3*v1*x4*y2 - u3*v2*x1*y4 + u3*v2*x4*y1 + u3*v4*x1*y2 - u3*v4*x2*y1 - u4*v1*x2*y3 + u4*v1*x3*y2 + u4*v2*x1*y3 - u4*v2*x3*y1 - u4*v3*x1*y2 + u4*v3*x2*y1)
    Cxy = (u1*u2p*v3*x4 - u1*u2p*v4*x3 - u1*u3p*v2*x4 + u1*u3p*v4*x2 + u1*u4p*v2*x3 - u1*u4p*v3*x2 - u1p*u2*v3*x4 + u1p*u2*v4*x3 + u1p*u3*v2*x4 - u1p*u3*v4*x2 - u1p*u4*v2*x3 + u1p*u4*v3*x2 + u2*u3p*v1*x4 - u2*u3p*v4*x1 - u2*u4p*v1*x3 + u2*u4p*v3*x1 - u2p*u3*v1*x4 + u2p*u3*v4*x1 + u2p*u4*v1*x3 - u2p*u4*v3*x1 + u3*u4p*v1*x2 - u3*u4p*v2*x1 - u3p*u4*v1*x2 + u3p*u4*v2*x1)/(u1*v2*x3*y4 - u1*v2*x4*y3 - u1*v3*x2*y4 + u1*v3*x4*y2 + u1*v4*x2*y3 - u1*v4*x3*y2 - u2*v1*x3*y4 + u2*v1*x4*y3 + u2*v3*x1*y4 - u2*v3*x4*y1 - u2*v4*x1*y3 + u2*v4*x3*y1 + u3*v1*x2*y4 - u3*v1*x4*y2 - u3*v2*x1*y4 + u3*v2*x4*y1 + u3*v4*x1*y2 - u3*v4*x2*y1 - u4*v1*x2*y3 + u4*v1*x3*y2 + u4*v2*x1*y3 - u4*v2*x3*y1 - u4*v3*x1*y2 + u4*v3*x2*y1)
    Dxx = (u1p*v2*x3*y4 - u1p*v2*x4*y3 - u1p*v3*x2*y4 + u1p*v3*x4*y2 + u1p*v4*x2*y3 - u1p*v4*x3*y2 - u2p*v1*x3*y4 + u2p*v1*x4*y3 + u2p*v3*x1*y4 - u2p*v3*x4*y1 - u2p*v4*x1*y3 + u2p*v4*x3*y1 + u3p*v1*x2*y4 - u3p*v1*x4*y2 - u3p*v2*x1*y4 + u3p*v2*x4*y1 + u3p*v4*x1*y2 - u3p*v4*x2*y1 - u4p*v1*x2*y3 + u4p*v1*x3*y2 + u4p*v2*x1*y3 - u4p*v2*x3*y1 - u4p*v3*x1*y2 + u4p*v3*x2*y1)/(u1*v2*x3*y4 - u1*v2*x4*y3 - u1*v3*x2*y4 + u1*v3*x4*y2 + u1*v4*x2*y3 - u1*v4*x3*y2 - u2*v1*x3*y4 + u2*v1*x4*y3 + u2*v3*x1*y4 - u2*v3*x4*y1 - u2*v4*x1*y3 + u2*v4*x3*y1 + u3*v1*x2*y4 - u3*v1*x4*y2 - u3*v2*x1*y4 + u3*v2*x4*y1 + u3*v4*x1*y2 - u3*v4*x2*y1 - u4*v1*x2*y3 + u4*v1*x3*y2 + u4*v2*x1*y3 - u4*v2*x3*y1 - u4*v3*x1*y2 + u4*v3*x2*y1)
    Dxy = (u1*u2p*x3*y4 - u1*u2p*x4*y3 - u1*u3p*x2*y4 + u1*u3p*x4*y2 + u1*u4p*x2*y3 - u1*u4p*x3*y2 - u1p*u2*x3*y4 + u1p*u2*x4*y3 + u1p*u3*x2*y4 - u1p*u3*x4*y2 - u1p*u4*x2*y3 + u1p*u4*x3*y2 + u2*u3p*x1*y4 - u2*u3p*x4*y1 - u2*u4p*x1*y3 + u2*u4p*x3*y1 - u2p*u3*x1*y4 + u2p*u3*x4*y1 + u2p*u4*x1*y3 - u2p*u4*x3*y1 + u3*u4p*x1*y2 - u3*u4p*x2*y1 - u3p*u4*x1*y2 + u3p*u4*x2*y1)/(u1*v2*x3*y4 - u1*v2*x4*y3 - u1*v3*x2*y4 + u1*v3*x4*y2 + u1*v4*x2*y3 - u1*v4*x3*y2 - u2*v1*x3*y4 + u2*v1*x4*y3 + u2*v3*x1*y4 - u2*v3*x4*y1 - u2*v4*x1*y3 + u2*v4*x3*y1 + u3*v1*x2*y4 - u3*v1*x4*y2 - u3*v2*x1*y4 + u3*v2*x4*y1 + u3*v4*x1*y2 - u3*v4*x2*y1 - u4*v1*x2*y3 + u4*v1*x3*y2 + u4*v2*x1*y3 - u4*v2*x3*y1 - u4*v3*x1*y2 + u4*v3*x2*y1)
    
    Cyx = (u1*v2*v3p*y4 - u1*v2*v4p*y3 - u1*v2p*v3*y4 + u1*v2p*v4*y3 + u1*v3*v4p*y2 - u1*v3p*v4*y2 - u2*v1*v3p*y4 + u2*v1*v4p*y3 + u2*v1p*v3*y4 - u2*v1p*v4*y3 - u2*v3*v4p*y1 + u2*v3p*v4*y1 + u3*v1*v2p*y4 - u3*v1*v4p*y2 - u3*v1p*v2*y4 + u3*v1p*v4*y2 + u3*v2*v4p*y1 - u3*v2p*v4*y1 - u4*v1*v2p*y3 + u4*v1*v3p*y2 + u4*v1p*v2*y3 - u4*v1p*v3*y2 - u4*v2*v3p*y1 + u4*v2p*v3*y1)/(u1*v2*x3*y4 - u1*v2*x4*y3 - u1*v3*x2*y4 + u1*v3*x4*y2 + u1*v4*x2*y3 - u1*v4*x3*y2 - u2*v1*x3*y4 + u2*v1*x4*y3 + u2*v3*x1*y4 - u2*v3*x4*y1 - u2*v4*x1*y3 + u2*v4*x3*y1 + u3*v1*x2*y4 - u3*v1*x4*y2 - u3*v2*x1*y4 + u3*v2*x4*y1 + u3*v4*x1*y2 - u3*v4*x2*y1 - u4*v1*x2*y3 + u4*v1*x3*y2 + u4*v2*x1*y3 - u4*v2*x3*y1 - u4*v3*x1*y2 + u4*v3*x2*y1)
    Cyy = (-u1*v2*v3p*x4 + u1*v2*v4p*x3 + u1*v2p*v3*x4 - u1*v2p*v4*x3 - u1*v3*v4p*x2 + u1*v3p*v4*x2 + u2*v1*v3p*x4 - u2*v1*v4p*x3 - u2*v1p*v3*x4 + u2*v1p*v4*x3 + u2*v3*v4p*x1 - u2*v3p*v4*x1 - u3*v1*v2p*x4 + u3*v1*v4p*x2 + u3*v1p*v2*x4 - u3*v1p*v4*x2 - u3*v2*v4p*x1 + u3*v2p*v4*x1 + u4*v1*v2p*x3 - u4*v1*v3p*x2 - u4*v1p*v2*x3 + u4*v1p*v3*x2 + u4*v2*v3p*x1 - u4*v2p*v3*x1)/(u1*v2*x3*y4 - u1*v2*x4*y3 - u1*v3*x2*y4 + u1*v3*x4*y2 + u1*v4*x2*y3 - u1*v4*x3*y2 - u2*v1*x3*y4 + u2*v1*x4*y3 + u2*v3*x1*y4 - u2*v3*x4*y1 - u2*v4*x1*y3 + u2*v4*x3*y1 + u3*v1*x2*y4 - u3*v1*x4*y2 - u3*v2*x1*y4 + u3*v2*x4*y1 + u3*v4*x1*y2 - u3*v4*x2*y1 - u4*v1*x2*y3 + u4*v1*x3*y2 + u4*v2*x1*y3 - u4*v2*x3*y1 - u4*v3*x1*y2 + u4*v3*x2*y1)
    Dyx = (-v1*v2p*x3*y4 + v1*v2p*x4*y3 + v1*v3p*x2*y4 - v1*v3p*x4*y2 - v1*v4p*x2*y3 + v1*v4p*x3*y2 + v1p*v2*x3*y4 - v1p*v2*x4*y3 - v1p*v3*x2*y4 + v1p*v3*x4*y2 + v1p*v4*x2*y3 - v1p*v4*x3*y2 - v2*v3p*x1*y4 + v2*v3p*x4*y1 + v2*v4p*x1*y3 - v2*v4p*x3*y1 + v2p*v3*x1*y4 - v2p*v3*x4*y1 - v2p*v4*x1*y3 + v2p*v4*x3*y1 - v3*v4p*x1*y2 + v3*v4p*x2*y1 + v3p*v4*x1*y2 - v3p*v4*x2*y1)/(u1*v2*x3*y4 - u1*v2*x4*y3 - u1*v3*x2*y4 + u1*v3*x4*y2 + u1*v4*x2*y3 - u1*v4*x3*y2 - u2*v1*x3*y4 + u2*v1*x4*y3 + u2*v3*x1*y4 - u2*v3*x4*y1 - u2*v4*x1*y3 + u2*v4*x3*y1 + u3*v1*x2*y4 - u3*v1*x4*y2 - u3*v2*x1*y4 + u3*v2*x4*y1 + u3*v4*x1*y2 - u3*v4*x2*y1 - u4*v1*x2*y3 + u4*v1*x3*y2 + u4*v2*x1*y3 - u4*v2*x3*y1 - u4*v3*x1*y2 + u4*v3*x2*y1)
    Dyy = (u1*v2p*x3*y4 - u1*v2p*x4*y3 - u1*v3p*x2*y4 + u1*v3p*x4*y2 + u1*v4p*x2*y3 - u1*v4p*x3*y2 - u2*v1p*x3*y4 + u2*v1p*x4*y3 + u2*v3p*x1*y4 - u2*v3p*x4*y1 - u2*v4p*x1*y3 + u2*v4p*x3*y1 + u3*v1p*x2*y4 - u3*v1p*x4*y2 - u3*v2p*x1*y4 + u3*v2p*x4*y1 + u3*v4p*x1*y2 - u3*v4p*x2*y1 - u4*v1p*x2*y3 + u4*v1p*x3*y2 + u4*v2p*x1*y3 - u4*v2p*x3*y1 - u4*v3p*x1*y2 + u4*v3p*x2*y1)/(u1*v2*x3*y4 - u1*v2*x4*y3 - u1*v3*x2*y4 + u1*v3*x4*y2 + u1*v4*x2*y3 - u1*v4*x3*y2 - u2*v1*x3*y4 + u2*v1*x4*y3 + u2*v3*x1*y4 - u2*v3*x4*y1 - u2*v4*x1*y3 + u2*v4*x3*y1 + u3*v1*x2*y4 - u3*v1*x4*y2 - u3*v2*x1*y4 + u3*v2*x4*y1 + u3*v4*x1*y2 - u3*v4*x2*y1 - u4*v1*x2*y3 + u4*v1*x3*y2 + u4*v2*x1*y3 - u4*v2*x3*y1 - u4*v3*x1*y2 + u4*v3*x2*y1)
    
    return np.array([[Axx,Axy,Bxx,Bxy],
                     [Ayx,Ayy,Byx,Byy],
                     [Cxx,Cxy,Dxx,Dxy],
                     [Cyx,Cyy,Dyx,Dyy]])

# Needs 2 raysets, before and after propagation - should call the orthogonal function
def compute_offdiagonal_abcd(ray1in,ray1out,ray2in,ray2out):
    
    x1 = ray1in[0,:]
    xp1 = ray1out[0,:]
    thx1 = ray1in[2,:]
    thxp1 = ray1in[2,:]
    
    y1 = ray1in[1,:]
    yp1 = ray1out[1,:]
    thy1 = ray1in[3,:]
    thyp1 = ray1in[3,:]
    
    x2 = ray2in[0,:]
    xp2 = ray2out[0,:]
    thx2 = ray2in[2,:]
    thxp2 = ray2in[2,:]
    
    y2 = ray2in[1,:]
    yp2 = ray2out[1,:]
    thy2 = ray2in[3,:]
    thyp2 = ray2in[3,:]

    # Compute the on-diagonals
    Axx,Bxx,Cxx,Dxx,Ayy,Byy,Cyy,Dyy = compute_orthogonal_abcd(ray1in,ray1out,ray2in,ray2out)

    # Compute the off-diagonals
    
    # A params
    Axy = (thy2*xp1 + Axx*thy1*x2 + Bxx*thx2*thy1 - thy1*xp2 - Bxx*thx1*thy2 - Axx*thy2*x1) / (-(thy1*y2 - thy2*y1))
    Ayx = (thx2*yp1 + Byy*thx1*thy2 + Ayy*thx1*y2 - thx1*yp2 - Byy*thx2*thy1 - Ayy*thx2*y1) / (-(thx1*x2 - thx2*x1))
    
    # B params
    Bxy = (xp1*y2 + Bxx*thx2*y1 + Axx*x2*y1 - xp2*y1 - Axx*x1*y2 - Bxx*thx1*y2) / (thy1*y2 - thy2*y1)
    Byx = (x2*yp1 + Ayy*x1*y2 + Byy*thy2*x1 - x1*yp2 - Byy*thy1*x2 - Ayy*x2*y1) / (thx1*x2 - thx2*x1)
    
    # C params
    Cxy = (thxp1*thy2 + Cxx*thy1*x2 + Dxx*thx2*thy1 - thxp2*thy1 - Cxx*thy2*x1 - Dxx*thx1*thy2) / (-(thy1*y2 - thy2*y1))
    Cyx = (thx2*thyp1 + Dyy*thx1*thy2 + Cyy*thx1*y2 - thx1*thyp2 - Cyy*thx2*y1 - Dyy*thx2*thy1) / (-(thx1*x2 - thx2*x1))
    
    # D params
    Dxy = (thxp1*y2 + Dxx*thx2*y1 + Cxx*x2*y1 - thxp2*y1 - Cxx*x1*y2 - Dxx*thx1*y2) / (thy1*y2 - thy2*y1)
    Dyx = (thyp1*x2 + Dyy*thy2*x1 + Cyy*x1*y2 - thyp2*x1 - Cyy*x2*y1 - Dyy*thy1*x2) / (thx1*x2 - thx2*x1)

    return np.array([[Axx,Axy,Bxx,Bxy],
                     [Ayx,Ayy,Byx,Byy],
                     [Cxx,Cxy,Dxx,Dxy],
                     [Cyx,Cyy,Dyx,Dyy]])

"""
Engine for interpretation of ray data as ray transfer matrices
"""

# formatted_rays = [x,y,thx,thy],[xp,yp,thxp,thyp]

# Interpreter for given ray format, begin with OpticStudio
def readrays(data):
    return formatted_rays

# Needs 2 raysets, before and after propagation
def compute_orthogonal_abcd(formatted_rays):

    A = ((r1p/th1) - (r2p/th2)) / ((r1/th1) - (r2/th2))
    B = ((r1p/r1) - (r2p/r2)) / ((th1/r1) - (th2/r2))
    C = ((th1p/th1) - (th2p/th2)) / ((r1/th1) - (r2/th2))
    D = ((th1p/r1) - (th2p/r2)) / ((th1/r1) - (th2/r2))

    return A,B,C,D

# Needs 2 raysets, before and after propagation - should call the orthogonal function
def compute_offdiagonal_abcd(formatted_rays):

    # Compute the on-diagonals
    Axx,Bxx,Cxx,Dxx = compute_orthogonal_abcd(formatted_rays)
    Ayy,Byy,Cyy,Dyy = compute_orthogonal_abcd(formatted_rays)

    # Compute the off-diagonals
    Axy = (thy2*xp1 + Axx*thy1*x2 + Bxx*thx2*thy1 - thy1*xp2 - Bxx*thx1*thy2 - Axx*thy2*x1) / (-(thy1*y2 - thy2*y1))
    Ayx = (thx2*yp1 + Byy*thx1*thy2 + Ayy*thx1*y2 - thx1*yp2 - Byy*thx2*thy1 - Ayy*thx2*y1) / (-(thx1*x2 - thx2*x1))
    Bxy = (xp1*y2 + Bxx*thx2*y1 + Axx*x2*y1 - xp2*y1 - Axx*x1*y2 - Bxx*thx1*y2) / (thy1*y2 - thy2*y1)
    Byx = (x2*yp1 + Ayy*x1*y2 + Byy*thy2*x1 - x1*yp2 - Byy*thy1*x2 - Ayy*x2*y1) / (thx1*x2 - thx2*x1)
    Cxy = (thxp1*thy2 + Cxx*thy1*x2 + Dxx*thx2*thy1 - thxp2*thy1 - Cxx*thy2*x1 - Dxx*thx1*thy2) / (-(thy1*y2 - thy2*y1))
    Cyx = (thx2*thyp1 + Dyy*thx1*thy2 + Cyy*thx1*y2 - thx1*thyp2 - Cyy*thx2*y1 - Dyy*thx2*thy1) / (-(thx1*x2 - thx2*x1))
    Dxy = (thxp1*y2 + Dxx*thx2*y1 + Cxx*x2*y1 - thxp2*y1 - Cxx*x1*y2 - Dxx*thx1*y2) / (thy1*y2 - thy2*y1)
    Dyx = (thyp1*x2 + Dyy*thy2*x1 + Cyy*x1*y2 - thyp2*x1 - Cyy*x2*y1 - Dyy*thy1*x2) / (thx1*x2 - thx2*x1)

    return np.array([[Axx,Axy,Bxx,Bxy],
                     [Ayx,Ayy,Byx,Byy],
                     [Cxx,Cxy,Dxx,Dyy],
                     [Cyx,Cyy,Dyx,Dyy]])

"""
Engine for interpretation of ray data as ray transfer matrices
"""

# formatted_rays = [x,y,thx,thy],[xp,yp,thxp,thyp]

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







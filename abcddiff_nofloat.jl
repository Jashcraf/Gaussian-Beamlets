using Symbolics

@variables x1 y1 thx1 thy1 xp1 yp1 thxp1 thyp1
@variables x2 y2 thx2 thy2 xp2 yp2 thxp2 thyp2
@variables x3 y3 thx3 thy3 xp3 yp3 thxp3 thyp3
@variables x4 y4 thx4 thy4 xp4 yp4 thxp4 thyp4
@variables Axx Axy Ayx Ayy Bxx Bxy Byx Byy Cxx Cxy Cyx Cyy Dxx Dxy Dyx Dyy

# Differential Ray Trace Data
aAB1 = Axx*x1 + Axy*y1 + Bxx*thx1 + Bxy*thy1 ~ xp1
bAB1 = Ayx*x1 + Ayy*y1 + Byx*thx1 + Byy*thy1 ~ yp1
aCD1 = Cxx*x1 + Cxy*y1 + Dxx*thx1 + Dxy*thy1 ~ thxp1
bCD1 = Cyx*x1 + Cyy*y1 + Dyx*thx1 + Dyy*thy1 ~ thyp1

aAB2 = Axx*x2 + Axy*y2 + Bxx*thx2 + Bxy*thy2 ~ xp2
bAB2 = Ayx*x2 + Ayy*y2 + Byx*thx2 + Byy*thy2 ~ yp2
aCD2 = Cxx*x2 + Cxy*y2 + Dxx*thx2 + Dxy*thy2 ~ thxp2
bCD2 = Cyx*x2 + Cyy*y2 + Dyx*thx2 + Dyy*thy2 ~ thyp2

aAB3 = Axx*x3 + Axy*y3 + Bxx*thx3 + Bxy*thy3 ~ xp3
bAB3 = Ayx*x3 + Ayy*y3 + Byx*thx3 + Byy*thy3 ~ yp3
aCD3 = Cxx*x3 + Cxy*y3 + Dxx*thx3 + Dxy*thy3 ~ thxp3
bCD3 = Cyx*x3 + Cyy*y3 + Dyx*thx3 + Dyy*thy3 ~ thyp3

aAB4 = Axx*x4 + Axy*y4 + Bxx*thx4 + Bxy*thy4 ~ xp4
bAB4 = Ayx*x4 + Ayy*y4 + Byx*thx4 + Byy*thy4 ~ yp4
aCD4 = Cxx*x4 + Cxy*y4 + Dxx*thx4 + Dxy*thy4 ~ thxp4
bCD4 = Cyx*x4 + Cyy*y4 + Dyx*thx4 + Dyy*thy4 ~ thyp4

#Ax = ((xp1/thx1)-(xp2/thy2))/((x1/thy1)-(x2/thx2)) ~ Axx
#Ay = ((yp1/thy1)-(yp2/thy2))/((y1/thy1)-(x2/thy2)) ~ Ayy

#Bx = ((xp1/x1)-(xp2/x2))/((thx1/x1)-(thx2/x2)) ~ Bxx
#By = ((yp1/y1)-(yp2/y2))/((thy1/y1)-(thy2/y2)) ~ Byy

#Cx = ((thxp1/thx1) - (thxp2/thx2))/((x1/thx1)-(x2/thx2)) ~ Cxx
#Cy = ((thyp1/thy1) - (thyp2/thy2))/((y1/thy1)-(y2/thy2)) ~ Cyy

#Dx = ((thxp1/x1)-(thxp2/x2))/((thx1/x1)-(thx2/x2)) ~ Dxx
#Dy = ((thyp1/y1)-(thyp2/y2))/((thy1/y1)-(thy2/y2)) ~ Dyy

# From Symplectic Condition
#det0 = Axx*Dxx + Axy*Dxy ~ Bxx*Cxx + Bxy*Cxy
#det1 = Axx*Dyx + Axy*Dyy ~ Bxx*Cyx + Bxy*Cyy
#det2 = Ayx*Dxx + Ayy*Dxy ~ Byx*Cxx + Byy*Cxy
#det3 = Ayx*Dyx + Ayy*Dyy ~ Byx*Cyx + Byy*Cyy

## I think these are redundant
# det4 = -Bxx*Cxx - Byx*Cyx + Axx*Dxx + Ayx*Dyx ~ sym(1)
# det5 = -Bxx*Cxy - Byx*Cyy + Axx*Dxy + Ayx*Dyy ~ sym(0)
# det6 = -Bxy*Cxx - Byy*Cyx + Axy*Dxx + Ayy*Dyx ~ sym(0)
# det7 = -Bxy*Cxy - Byy*Cyy + Axy*Dxy + Ayy*Dyy ~ sym(1)

#cAB = Axx*Byx + Axy*Byy ~ Ayx*Bxx + Ayy*Bxy
#aBD = Bxx*Dxy + Byx*Dyy ~ Byy*Dyx + Bxy*Dxx
#cCD = Cyx*Dxx + Cyy*Dxy ~ Cxy*Dyy + Cxx*Dyx
#aAC = Ayx*Cxx + Ayy*Cyx ~ Ayx*Cyy + Axx*Cxy
# dAB = -Axx*Byx - Axy*Byy + Ayx*Bxx + Ayy*Bxy ~ sym(0)
# bBD = -Bxx*Dxy + Bxy*Dxx - Byx*Dyy + Byy*Dyx ~ sym(0)
# dCD = Cxx*Dyx + Cxy*Dyy - Cyx*Dxx - Cyy*Dxy ~ sym(0)
# bAC = Axx*Cxy - Ayx*Cxx + Ayx*Cyy - Ayy*Cyx ~ sym(0)
#abcd1 = Axx*Dxx - Bxx*Cxx ~ 1
#abcd2 = Axx*x1 + Bxx*thx1 ~ xp1
#abcd3 = Cxx*x1 + Dxx*thx1 ~ thxp1
#abcd4 = Axx*x2 + Bxx*thx2 ~ xp2
# abcd5 = Cxx*x2 + Dxx*thx2 ~ thxp2

# Symbolics.solve_for([abcd1,abcd2,abcd3,abcd4],[Axx,Bxx,Cxx,Dxx])
# aAB2,bAB2,aCD2,bCD2,
# det0,det1,det2,det3,cAB,aBD,cCD,aAC
#Symbolics.solve_for([aAB1,bAB1,aCD1,bCD1,aAB2,bAB2,aCD2,bCD2,aAB3,bAB3,aCD3,bCD3,aAB4,bAB4,aCD4,bCD4],
#                    [Axx,Axy,Ayx,Ayy,Bxx,Bxy,Byx,Byy,Cxx,Cxy,Cyx,Cyy,Dxx,Dxy,Dyx,Dyy],simplify=true)

# Solve for A and B first row
# Symbolics.solve_for([aAB1,aAB2,aAB3,aAB4],
# 		   [Axx,Axy,Bxx,Bxy])

# Solve for A and B second row
#Symbolics.solve_for([bAB1,bAB2,bAB3,bAB4],
# 		   [Ayx,Ayy,Byx,Byy])

# Solve for C and D first row
#Symbolics.solve_for([aCD1,aCD2,aCD3,aCD4],
# 		   [Cxx,Cxy,Dxx,Dxy])

# Solve for C and D second row
Symbolics.solve_for([bCD1,bCD2,bCD3,bCD4],
 		   [Cyx,Cyy,Dyx,Dyy])
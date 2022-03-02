#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 15:37:35 2021

@author: jashcraft
"""

import numpy as np
import differential as dif

# px = np.linspace(-1200,1200,50+1)
# px,py = np.meshgrid(px,px)
# px = np.ravel(px)
# py = np.ravel(py)

# deltaray1x = 1e-3
# deltaang1x = np.pi/2 - np.arccos(deltaray1x)
# deltaray1y = 2e-3
# deltaang1y = np.pi/2 - np.arccos(deltaray1y)

# deltaray2x = -3e-3
# deltaang2x = np.pi/2 - np.arccos(deltaray2x)
# deltaray2y = 4e-3
# deltaang2y = np.pi/2 - np.arccos(deltaray2y)

# deltaray3x = -5e-3
# deltaang3x = np.pi/2 - np.arccos(deltaray3x)
# deltaray3y = -6e-3
# deltaang3y = np.pi/2 - np.arccos(deltaray3y)

# deltaray4x = 7e-3
# deltaang4x = np.pi/2 - np.arccos(deltaray4x)
# deltaray4y = -8e-3
# deltaang4y = np.pi/2 - np.arccos(deltaray4y)

# raydelta1 = 0.00005
# raydeltang1 = np.pi/2 - np.arccos(raydelta)
# raydelta1 = 0.00005
# raydeltang1 = np.pi/2 - np.arccos(raydelta)
# raydelta1 = 0.00005
# raydeltang1 = np.pi/2 - np.arccos(raydelta)
# raydelta1 = 0.00005
# raydeltang1 = np.pi/2 - np.arccos(raydelta)

# ray0in = np.array([px,
#                    py,
#                    0*px,
#                    0*py])

# ray1in = np.array([px+deltaray1x,
#                    py+deltaray1y,
#                    0*px+deltaang1x,
#                    0*py+deltaang1y])

# ray2in = np.array([px+deltaray2x,
#                    py+deltaray2y,
#                    0*px+deltaang2x,
#                    0*py+deltaang2y])

# ray3in = np.array([px+deltaray3x,
#                    py+deltaray3y,
#                    0*px+deltaang3x,
#                    0*py+deltaang3y])

# ray4in = np.array([px+deltaray4x,
#                    py+deltaray4y,
#                    0*px+deltaang4x,
#                    0*py+deltaang4y])

# rayout1x = 'parabola-test/parabolanom.zmx_50_5e-05_5e-05_3_xray.txt'
# rayout1y = 'parabola-test/parabolanom.zmx_50_5e-05_5e-05_3_yray.txt'
# rayout1a = 'parabola-test/parabolanom.zmx_50_5e-05_5e-05_3_aray.txt'
# rayout1b = 'parabola-test/parabolanom.zmx_50_5e-05_5e-05_3_bray.txt'

# rayout2x = 'parabola-test/parabolanom.zmx_50_-5e-05_5e-05_3_xray.txt'
# rayout2y = 'parabola-test/parabolanom.zmx_50_-5e-05_5e-05_3_yray.txt'
# rayout2a = 'parabola-test/parabolanom.zmx_50_-5e-05_5e-05_3_aray.txt'
# rayout2b = 'parabola-test/parabolanom.zmx_50_-5e-05_5e-05_3_bray.txt'

# rayout3x = 'parabola-test/parabolanom.zmx_50_5e-05_-5e-05_3_xray.txt'
# rayout3y = 'parabola-test/parabolanom.zmx_50_5e-05_-5e-05_3_yray.txt'
# rayout3a = 'parabola-test/parabolanom.zmx_50_5e-05_-5e-05_3_aray.txt'
# rayout3b = 'parabola-test/parabolanom.zmx_50_5e-05_-5e-05_3_bray.txt'

# rayout4x = 'parabola-test/parabolanom.zmx_50_-5e-05_-5e-05_3_xray.txt'
# rayout4y = 'parabola-test/parabolanom.zmx_50_-5e-05_-5e-05_3_yray.txt'
# rayout4a = 'parabola-test/parabolanom.zmx_50_-5e-05_-5e-05_3_aray.txt'
# rayout4b = 'parabola-test/parabolanom.zmx_50_-5e-05_-5e-05_3_bray.txt'

# -- - --- 

# rayout1x = 'parabola-test/parabolanom.zmx_50_1e-05_2e-05_3_xray.txt'
# rayout1y = 'parabola-test/parabolanom.zmx_50_1e-05_2e-05_3_yray.txt'
# rayout1a = 'parabola-test/parabolanom.zmx_50_1e-05_2e-05_3_aray.txt'
# rayout1b = 'parabola-test/parabolanom.zmx_50_1e-05_2e-05_3_bray.txt'

# rayout2x = 'parabola-test/parabolanom.zmx_50_-3e-05_4e-05_3_xray.txt'
# rayout2y = 'parabola-test/parabolanom.zmx_50_-3e-05_4e-05_3_yray.txt'
# rayout2a = 'parabola-test/parabolanom.zmx_50_-3e-05_4e-05_3_aray.txt'
# rayout2b = 'parabola-test/parabolanom.zmx_50_-3e-05_4e-05_3_bray.txt'

# rayout3x = 'parabola-test/parabolanom.zmx_50_-5e-05_-6e-05_3_xray.txt'
# rayout3y = 'parabola-test/parabolanom.zmx_50_-5e-05_-6e-05_3_yray.txt'
# rayout3a = 'parabola-test/parabolanom.zmx_50_-5e-05_-6e-05_3_aray.txt'
# rayout3b = 'parabola-test/parabolanom.zmx_50_-5e-05_-6e-05_3_bray.txt'

# rayout4x = 'parabola-test/parabolanom.zmx_50_7e-05_-8e-05_3_xray.txt'
# rayout4y = 'parabola-test/parabolanom.zmx_50_7e-05_-8e-05_3_yray.txt'
# rayout4a = 'parabola-test/parabolanom.zmx_50_7e-05_-8e-05_3_aray.txt'
# rayout4b = 'parabola-test/parabolanom.zmx_50_7e-05_-8e-05_3_bray.txt'

# rayout1x = 'parabola-test/freespace.zmx_50_0.001_0.002_3_xray.txt'
# rayout1y = 'parabola-test/freespace.zmx_50_0.001_0.002_3_yray.txt'
# rayout1a = 'parabola-test/freespace.zmx_50_0.001_0.002_3_aray.txt'
# rayout1b = 'parabola-test/freespace.zmx_50_0.001_0.002_3_bray.txt'

# rayout2x = 'parabola-test/freespace.zmx_50_-0.003_0.004_3_xray.txt'
# rayout2y = 'parabola-test/freespace.zmx_50_-0.003_0.004_3_yray.txt'
# rayout2a = 'parabola-test/freespace.zmx_50_-0.003_0.004_3_aray.txt'
# rayout2b = 'parabola-test/freespace.zmx_50_-0.003_0.004_3_bray.txt'

# rayout3x = 'parabola-test/freespace.zmx_50_-0.005_-0.006_3_xray.txt'
# rayout3y = 'parabola-test/freespace.zmx_50_-0.005_-0.006_3_yray.txt'
# rayout3a = 'parabola-test/freespace.zmx_50_-0.005_-0.006_3_aray.txt'
# rayout3b = 'parabola-test/freespace.zmx_50_-0.005_-0.006_3_bray.txt'

# rayout4x = 'parabola-test/freespace.zmx_50_0.007_-0.008_3_xray.txt'
# rayout4y = 'parabola-test/freespace.zmx_50_0.007_-0.008_3_yray.txt'
# rayout4a = 'parabola-test/freespace.zmx_50_0.007_-0.008_3_aray.txt'
# rayout4b = 'parabola-test/freespace.zmx_50_0.007_-0.008_3_bray.txt'

# # ----
# rayin0x  = 'parabola-nom/parabolanom.zmx_50_0_0_0_xray.txt'
# rayin0y  = 'parabola-nom/parabolanom.zmx_50_0_0_0_yray.txt'
# rayin0a  = 'parabola-nom/parabolanom.zmx_50_0_0_0_aray.txt'
# rayin0b  = 'parabola-nom/parabolanom.zmx_50_0_0_0_bray.txt'

# rayout0x = 'parabola-nom/parabolanom.zmx_50_0_0_3_xray.txt'
# rayout0y = 'parabola-nom/parabolanom.zmx_50_0_0_3_yray.txt'
# rayout0a = 'parabola-nom/parabolanom.zmx_50_0_0_3_aray.txt'
# rayout0b = 'parabola-nom/parabolanom.zmx_50_0_0_3_bray.txt'

# rayin1x  = 'parabola-nom/parabolanom.zmx_50_0.001_0.002_0_xray.txt'
# rayin1y  = 'parabola-nom/parabolanom.zmx_50_0.001_0.002_0_yray.txt'
# rayin1a  = 'parabola-nom/parabolanom.zmx_50_0.001_0.002_0_aray.txt'
# rayin1b  = 'parabola-nom/parabolanom.zmx_50_0.001_0.002_0_bray.txt'

# rayout1x = 'parabola-nom/parabolanom.zmx_50_0.001_0.002_3_xray.txt'
# rayout1y = 'parabola-nom/parabolanom.zmx_50_0.001_0.002_3_yray.txt'
# rayout1a = 'parabola-nom/parabolanom.zmx_50_0.001_0.002_3_aray.txt'
# rayout1b = 'parabola-nom/parabolanom.zmx_50_0.001_0.002_3_bray.txt'

# rayin2x  = 'parabola-nom/parabolanom.zmx_50_-0.003_0.004_0_xray.txt'
# rayin2y  = 'parabola-nom/parabolanom.zmx_50_-0.003_0.004_0_yray.txt'
# rayin2a  = 'parabola-nom/parabolanom.zmx_50_-0.003_0.004_0_aray.txt'
# rayin2b  = 'parabola-nom/parabolanom.zmx_50_-0.003_0.004_0_bray.txt'

# rayout2x = 'parabola-nom/parabolanom.zmx_50_-0.003_0.004_3_xray.txt'
# rayout2y = 'parabola-nom/parabolanom.zmx_50_-0.003_0.004_3_yray.txt'
# rayout2a = 'parabola-nom/parabolanom.zmx_50_-0.003_0.004_3_aray.txt'
# rayout2b = 'parabola-nom/parabolanom.zmx_50_-0.003_0.004_3_bray.txt'

# rayin3x  = 'parabola-nom/parabolanom.zmx_50_-0.005_-0.006_0_xray.txt'
# rayin3y  = 'parabola-nom/parabolanom.zmx_50_-0.005_-0.006_0_yray.txt'
# rayin3a  = 'parabola-nom/parabolanom.zmx_50_-0.005_-0.006_0_aray.txt'
# rayin3b  = 'parabola-nom/parabolanom.zmx_50_-0.005_-0.006_0_bray.txt'

# rayout3x = 'parabola-nom/parabolanom.zmx_50_-0.005_-0.006_3_xray.txt'
# rayout3y = 'parabola-nom/parabolanom.zmx_50_-0.005_-0.006_3_yray.txt'
# rayout3a = 'parabola-nom/parabolanom.zmx_50_-0.005_-0.006_3_aray.txt'
# rayout3b = 'parabola-nom/parabolanom.zmx_50_-0.005_-0.006_3_bray.txt'

# rayin4x  = 'parabola-nom/parabolanom.zmx_50_0.007_-0.008_0_xray.txt'
# rayin4y  = 'parabola-nom/parabolanom.zmx_50_0.007_-0.008_0_yray.txt'
# rayin4a  = 'parabola-nom/parabolanom.zmx_50_0.007_-0.008_0_aray.txt'
# rayin4b  = 'parabola-nom/parabolanom.zmx_50_0.007_-0.008_0_bray.txt'

# rayout4x = 'parabola-nom/parabolanom.zmx_50_0.007_-0.008_3_xray.txt'
# rayout4y = 'parabola-nom/parabolanom.zmx_50_0.007_-0.008_3_yray.txt'
# rayout4a = 'parabola-nom/parabolanom.zmx_50_0.007_-0.008_3_aray.txt'
# rayout4b = 'parabola-nom/parabolanom.zmx_50_0.007_-0.008_3_bray.txt'

# ray1in = dif.readrays([rayin1x,rayin1y,rayin1a,rayin1b])
# ray2in = dif.readrays([rayin2x,rayin2y,rayin2a,rayin2b])

# ----
rayin0x  = 'direct-ray-test/paraxialnom.zmx_50_0_0_1_xray.txt'
rayin0y  = 'direct-ray-test/paraxialnom.zmx_50_0_0_1_yray.txt'
rayin0z  = 'direct-ray-test/paraxialnom.zmx_50_0_0_1_zray.txt'
rayin0a  = 'direct-ray-test/paraxialnom.zmx_50_0_0_1_aray.txt'
rayin0b  = 'direct-ray-test/paraxialnom.zmx_50_0_0_1_bray.txt'
rayin0g  = 'direct-ray-test/paraxialnom.zmx_50_0_0_1_gray.txt'

rayout0x = 'direct-ray-test/paraxialnom.zmx_50_0_0_3_xray.txt'
rayout0y = 'direct-ray-test/paraxialnom.zmx_50_0_0_3_yray.txt'
rayout0z = 'direct-ray-test/paraxialnom.zmx_50_0_0_3_zray.txt'
rayout0a = 'direct-ray-test/paraxialnom.zmx_50_0_0_3_aray.txt'
rayout0b = 'direct-ray-test/paraxialnom.zmx_50_0_0_3_bray.txt'
rayout0g = 'direct-ray-test/paraxialnom.zmx_50_0_0_3_gray.txt'

rayin1x  = 'direct-ray-test/paraxialnom.zmx_50_1e-06_2e-06_1_xray.txt'
rayin1y  = 'direct-ray-test/paraxialnom.zmx_50_1e-06_2e-06_1_yray.txt'
rayin1z  = 'direct-ray-test/paraxialnom.zmx_50_1e-06_2e-06_1_zray.txt'
rayin1a  = 'direct-ray-test/paraxialnom.zmx_50_1e-06_2e-06_1_aray.txt'
rayin1b  = 'direct-ray-test/paraxialnom.zmx_50_1e-06_2e-06_1_bray.txt'
rayin1g  = 'direct-ray-test/paraxialnom.zmx_50_1e-06_2e-06_1_gray.txt'

rayout1x = 'direct-ray-test/paraxialnom.zmx_50_1e-06_2e-06_3_xray.txt'
rayout1y = 'direct-ray-test/paraxialnom.zmx_50_1e-06_2e-06_3_yray.txt'
rayout1z = 'direct-ray-test/paraxialnom.zmx_50_1e-06_2e-06_3_zray.txt'
rayout1a = 'direct-ray-test/paraxialnom.zmx_50_1e-06_2e-06_3_aray.txt'
rayout1b = 'direct-ray-test/paraxialnom.zmx_50_1e-06_2e-06_3_bray.txt'
rayout1g = 'direct-ray-test/paraxialnom.zmx_50_1e-06_2e-06_3_gray.txt'

rayin2x  = 'direct-ray-test/paraxialnom.zmx_50_-3e-06_4e-06_1_xray.txt'
rayin2y  = 'direct-ray-test/paraxialnom.zmx_50_-3e-06_4e-06_1_yray.txt'
rayin2z  = 'direct-ray-test/paraxialnom.zmx_50_-3e-06_4e-06_1_zray.txt'
rayin2a  = 'direct-ray-test/paraxialnom.zmx_50_-3e-06_4e-06_1_aray.txt'
rayin2b  = 'direct-ray-test/paraxialnom.zmx_50_-3e-06_4e-06_1_bray.txt'
rayin2g  = 'direct-ray-test/paraxialnom.zmx_50_-3e-06_4e-06_1_gray.txt'

rayout2x = 'direct-ray-test/paraxialnom.zmx_50_-3e-06_4e-06_3_xray.txt'
rayout2y = 'direct-ray-test/paraxialnom.zmx_50_-3e-06_4e-06_3_yray.txt'
rayout2z = 'direct-ray-test/paraxialnom.zmx_50_-3e-06_4e-06_3_zray.txt'
rayout2a = 'direct-ray-test/paraxialnom.zmx_50_-3e-06_4e-06_3_aray.txt'
rayout2b = 'direct-ray-test/paraxialnom.zmx_50_-3e-06_4e-06_3_bray.txt'
rayout2g = 'direct-ray-test/paraxialnom.zmx_50_-3e-06_4e-06_3_gray.txt'

rayin3x  = 'direct-ray-test/paraxialnom.zmx_50_-5e-06_-6e-06_1_xray.txt'
rayin3y  = 'direct-ray-test/paraxialnom.zmx_50_-5e-06_-6e-06_1_yray.txt'
rayin3z  = 'direct-ray-test/paraxialnom.zmx_50_-5e-06_-6e-06_1_zray.txt'
rayin3a  = 'direct-ray-test/paraxialnom.zmx_50_-5e-06_-6e-06_1_aray.txt'
rayin3b  = 'direct-ray-test/paraxialnom.zmx_50_-5e-06_-6e-06_1_bray.txt'
rayin3g  = 'direct-ray-test/paraxialnom.zmx_50_-5e-06_-6e-06_1_gray.txt'

rayout3x = 'direct-ray-test/paraxialnom.zmx_50_-5e-06_-6e-06_3_xray.txt'
rayout3y = 'direct-ray-test/paraxialnom.zmx_50_-5e-06_-6e-06_3_yray.txt'
rayout3z = 'direct-ray-test/paraxialnom.zmx_50_-5e-06_-6e-06_3_zray.txt'
rayout3a = 'direct-ray-test/paraxialnom.zmx_50_-5e-06_-6e-06_3_aray.txt'
rayout3b = 'direct-ray-test/paraxialnom.zmx_50_-5e-06_-6e-06_3_bray.txt'
rayout3g = 'direct-ray-test/paraxialnom.zmx_50_-5e-06_-6e-06_3_gray.txt'

rayin4x  = 'direct-ray-test/paraxialnom.zmx_50_7e-06_-8e-06_1_xray.txt'
rayin4y  = 'direct-ray-test/paraxialnom.zmx_50_7e-06_-8e-06_1_yray.txt'
rayin4z  = 'direct-ray-test/paraxialnom.zmx_50_7e-06_-8e-06_1_zray.txt'
rayin4a  = 'direct-ray-test/paraxialnom.zmx_50_7e-06_-8e-06_1_aray.txt'
rayin4b  = 'direct-ray-test/paraxialnom.zmx_50_7e-06_-8e-06_1_bray.txt'
rayin4g  = 'direct-ray-test/paraxialnom.zmx_50_7e-06_-8e-06_1_gray.txt'

rayout4x = 'direct-ray-test/paraxialnom.zmx_50_7e-06_-8e-06_3_xray.txt'
rayout4y = 'direct-ray-test/paraxialnom.zmx_50_7e-06_-8e-06_3_yray.txt'
rayout4z = 'direct-ray-test/paraxialnom.zmx_50_7e-06_-8e-06_3_zray.txt'
rayout4a = 'direct-ray-test/paraxialnom.zmx_50_7e-06_-8e-06_3_aray.txt'
rayout4b = 'direct-ray-test/paraxialnom.zmx_50_7e-06_-8e-06_3_bray.txt'
rayout4g = 'direct-ray-test/paraxialnom.zmx_50_7e-06_-8e-06_3_gray.txt'

ray0in = dif.readrays([rayin0x,rayin0y,rayin0z,rayin0a,rayin0b,rayin0g])[:,:2500]
ray1in = dif.readrays([rayin1x,rayin1y,rayin1z,rayin1a,rayin1b,rayin1g])[:,:2500]
ray2in = dif.readrays([rayin2x,rayin2y,rayin2z,rayin2a,rayin2b,rayin2g])[:,:2500]
ray3in = dif.readrays([rayin3x,rayin3y,rayin3z,rayin3a,rayin3b,rayin3g])[:,:2500]
ray4in = dif.readrays([rayin4x,rayin4y,rayin4z,rayin4a,rayin4b,rayin4g])[:,:2500]

ray0out = dif.readrays([rayout0x,rayout0y,rayout0z,rayout0a,rayout0b,rayout0g])[:,:2500]
ray1out = dif.readrays([rayout1x,rayout1y,rayout1z,rayout1a,rayout1b,rayout1g])[:,:2500]
ray2out = dif.readrays([rayout2x,rayout2y,rayout2z,rayout2a,rayout2b,rayout2g])[:,:2500]
ray3out = dif.readrays([rayout3x,rayout3y,rayout3z,rayout3a,rayout3b,rayout3g])[:,:2500]
ray4out = dif.readrays([rayout4x,rayout4y,rayout4z,rayout4a,rayout4b,rayout4g])[:,:2500]

# bypass read ray data

x = np.linspace(-1,1,256)
x,y = np.meshgrid(x,x)
x = np.ravel(x)
y = np.ravel(y)

# ray0in = np.array([x+0e-6,y+0e-6,
#                    0*x+0e-4,0*y+0e-4])

# ray1in = np.array([x+1.1e-6,y+1.1e-6,
#                    0*x+1.1e-4,0*y+1.1e-4])
# ray2in = np.array([x+1.2e-6,y-1.2e-6,
#                    0*x+1.2e-4,0*y-1.2e-4])
# ray3in = np.array([x-1.3e-6,y+1.3e-6,
#                    0*x-1.3e-4,0*y+1.3e-4])
# ray4in = np.array([x-1.4e-6,y-1.4e-6,
#                    0*x-1.4e-4,0*y-1.4e-4])

# import gfuncs as gf
# lens = gf.ThinLens(-10,ray0in.shape[1])
# ray1out = gf.matmultrays(lens,ray1in)
# ray2out = gf.matmultrays(lens,ray2in)
# ray3out = gf.matmultrays(lens,ray3in)
# ray4out = gf.matmultrays(lens,ray4in)



# ray2in = ray1in-ray2in
# ray2out = ray1out-ray2out

# abcd = dif.compute_full_offdiagonal_abcd(ray1in-ray0in,ray1out-ray0out,
#                                           ray2in-ray0in,ray2out-ray0out,
#                                           ray3in-ray0in,ray3out-ray0out,
#                                           ray4in-ray0in,ray4out-ray0out)

# abcd = dif.compute_full_offdiagonal_abcd(ray1in+raydiff,ray1out+raydiff,
#                                          ray2in+raydiff,ray2out+raydiff,
#                                          ray3in+raydiff,ray3out+raydiff,
#                                          ray4in+raydiff,ray4out+raydiff)

abcd = dif.compute_full_offdiagonal_abcd(ray1in, ray1out, 
                                          ray2in, ray2out, 
                                          ray3in, ray3out, 
                                          ray4in, ray4out)

#abcd = # dif.compute_offdiagonal_abcd(ray1in,ray1out,ray2in,ray2out)

def plotrays(rayset):
    
    import matplotlib.pyplot as plt
    plt.figure(figsize=[14,7])
    plt.subplot(121)
    plt.scatter(rayset[0,:],rayset[1,:])
    plt.title('position')
    plt.subplot(122)
    plt.scatter(rayset[2,:],rayset[3,:])
    plt.title('angle')
    plt.show()
    
def rawcompare(r1,r2,r3,r4):
    import matplotlib.pyplot as plt
    
    plt.figure()
    plt.title('X Data')
    plt.plot(r1[0,:],label='ray 1')
    plt.plot(r2[0,:],label='ray 2')
    plt.plot(r3[0,:],label='ray 3')
    plt.plot(r4[0,:],label='ray 4')
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.title('Y Data')
    plt.plot(r1[1,:],label='ray 1')
    plt.plot(r2[1,:],label='ray 2')
    plt.plot(r3[1,:],label='ray 3')
    plt.plot(r4[1,:],label='ray 4')
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.title('A Data')
    plt.plot(r1[2,:],label='ray 1')
    plt.plot(r2[2,:],label='ray 2')
    plt.plot(r3[2,:],label='ray 3')
    plt.plot(r4[2,:],label='ray 4')
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.title('B Data')
    plt.plot(r1[3,:],label='ray 1')
    plt.plot(r2[3,:],label='ray 2')
    plt.plot(r3[3,:],label='ray 3')
    plt.plot(r4[3,:],label='ray 4')
    plt.legend()
    plt.show()

def compare_abcdels(Px,Py):
    
    import matplotlib.pyplot as plt
    plt.figure(figsize=[10,10])
    plt.scatter(Px,Py)
    #plt.plot(Px,label='X data')
    #plt.plot(Py,label='Y data')
    #plt.legend()
    #plt.ylim([0,2.3*2.4e3])
    plt.show()
    
def checksysfornan(sys):
    
    print('Axx nan values',sys[0,0,:][np.isnan(sys[0,0,:])])
    print('Axy nan values',sys[0,1,:][np.isnan(sys[0,1,:])])
    print('Ayx nan values',sys[1,0,:][np.isnan(sys[1,0,:])])
    print('Ayy nan values',sys[1,1,:][np.isnan(sys[1,1,:])])
    
    print('Bxx nan values',sys[0,2,:][np.isnan(sys[0,2,:])])
    print('Bxy nan values',sys[0,3,:][np.isnan(sys[0,3,:])])
    print('Byx nan values',sys[1,2,:][np.isnan(sys[1,2,:])])
    print('Byy nan values',sys[1,3,:][np.isnan(sys[1,3,:])])
    
    print('Cxx nan values',sys[2,0,:][np.isnan(sys[2,0,:])])
    print('Cxy nan values',sys[2,1,:][np.isnan(sys[2,1,:])])
    print('Cyx nan values',sys[3,0,:][np.isnan(sys[3,0,:])])
    print('Cyy nan values',sys[3,1,:][np.isnan(sys[3,1,:])])
    
    print('Dxx nan values',sys[2,2,:][np.isnan(sys[2,2,:])])
    print('Dxy nan values',sys[2,3,:][np.isnan(sys[2,3,:])])
    print('Dyx nan values',sys[3,2,:][np.isnan(sys[3,2,:])])
    print('Dyy nan values',sys[3,3,:][np.isnan(sys[3,3,:])])
    
import utils as ju
import gfuncs as gf

ju.fourbyfour(abcd,2,ray0in)
checksysfornan(abcd)
#rayouttest = gf.matmultrays(abcd,ray0in)
#field = gf.eval_gausfield(ray0in,abcd,1e-6,1e-3,512,amps=None)






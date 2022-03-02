#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 13:18:04 2022

@author: jashcraft
"""

import raytrace as rt
import gbd
import utils as ut
import poppy
import numpy as np
import astropy.units as u

wlen = 2.2e-6
size = 25.4e-3

rays = rt.makerays(size, 90,circle=False)
rayz = rt.ZernikeWFE(90, size, rays, 4, wlen)
nrays = int((rays.shape[1]))
lens = rt.ThinLens(10*size,nrays)
dist = rt.FreeSpace(10*size, nrays)

syst = rt.matmultlist(dist,lens)
rayo = rt.matmultrays(syst,rayz)

efie = np.reshape(gbd.eval_gausfield(rays, syst, wlen, 1e-3, 256),[256,256])
ut.amp_phase_plot(efie,vmin=1e-1)

"""
Now the Poppy system
"""

osys = poppy.FresnelOpticalSystem(pupil_diameter=size*u.m, npix=256, beam_ratio=0.25)
l1 = poppy.QuadraticLens(10*size)
im = poppy.ScalarTransmission()#poppy.Detector(pixelscale=size*u.m/256,fov_pixels=256)
#zn = poppy.ZernikeWFE()

osys.add_optic(l1)
osys.add_optic(im,distance=10*size*u.m)
import matplotlib.pyplot as plt
psf = osys.calc_psf(wavelength=wlen)
plt.figure()
poppy.display_psf(psf,imagecrop=1e-3,cmap='plasma')
plt.set_cmap('plasma')
plt.show

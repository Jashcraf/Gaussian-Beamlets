#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 10:33:01 2022

@author: jashcraft
"""

import numpy as np

def readrays(fn):
    
    raydata = np.loadtxt(fn)
    
    x = raydata[0]
    y = raydata[1]
    a = raydata[2]
    b = raydata[3]
    
    v = np.tan(np.arccos(a)-np.pi)
    u = np.tan(np.arccos(b)-np.pi)
    
    rays = np.array([x,
                     y,
                     u,
                     v])
    
    return rays

def construct_rtm(rayin1,rayout1,
                  rayin2,rayout2,
                  rayin3,rayout3,
                  rayin4,rayout4):
    
    
    rayin = np.array([[rayin1],
                      [rayin2],
                      [rayin3],
                      [rayin4]])
    
    rayout = np.array([[rayout1],
                       [rayout2],
                       [rayout3],
                       [rayout4]])
    
    rtm = np.matmul(rayout,np.linalg.inv(rayin)) 
    
    return rtm
    
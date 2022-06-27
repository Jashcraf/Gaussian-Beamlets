#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 10:33:01 2022

@author: jashcraft
"""

import numpy as np

def readrays(fn):
    
    raydata = np.loadtxt(fn,skiprows=1,delimiter=',')

    # print(raydata.shape)
    
    x = raydata[:,0]
    y = raydata[:,1]
    z = raydata[:,2]
    a = raydata[:,3]
    b = raydata[:,4]
    g = raydata[:,5]

    # r = np.array([x,y,z])
    # k = np.array([a,b,g])
    # norm = np.linalg.norm(r,axis=0)
    # kx = np.array([a,g])
    # xz = kx/np.linalg.norm(kx,axis=0)
    # xz = xz[1,:]
    # ky = np.array([b,g])
    # yz = ky/np.linalg.norm(ky,axis=0)
    # yz = yz[1,:]
    # print('norm = ',norm)
    
    v = np.tan(np.pi/2-np.arccos(b))
    u = np.tan(np.pi/2-np.arccos(a))
    # u = np.tan(np.arccos(a))
    # v = np.tan(np.arccos(b))
    
    rays = np.array([x,
                     y,
                     u,
                     v])

    # print(rays.shape)

    # rays = rays[:,(x**2 + y**2 + u**2 + v**2)**(1/2) > 0]
    
    return rays

def ConstructDiffMat(rayin0,rayout0,
                    rayin1,rayout1,
                    rayin2,rayout2,
                    rayin3,rayout3,
                    rayin4,rayout4):

    # Compute Differences
    rayin1 -= rayin0
    rayin2 -= rayin0
    rayin3 -= rayin0
    rayin4 -= rayin0

    rayout1 -= rayout0
    rayout2 -= rayout0
    rayout3 -= rayout0
    rayout4 -= rayout0

    # grab parameters
    xin1,yin1,uin1,vin1 = _grab_ray_data(rayin1)
    xin2,yin2,uin2,vin2 = _grab_ray_data(rayin2)
    xin3,yin3,uin3,vin3 = _grab_ray_data(rayin3)
    xin4,yin4,uin4,vin4 = _grab_ray_data(rayin4)

    xout1,yout1,uout1,vout1 = _grab_ray_data(rayout1)
    xout2,yout2,uout2,vout2 = _grab_ray_data(rayout2)
    xout3,yout3,uout3,vout3 = _grab_ray_data(rayout3)
    xout4,yout4,uout4,vout4 = _grab_ray_data(rayout4)

    # You only need some of these
    # 1 - the ray mod in +x position, relates to A matrix
    # 2 - the ray mod in +y position, relates to A Matrix
    # 3 - the ray mod in +x angle
    # 4 - the ray mod in +y angle

    # convert the parameters that matter
    # dx1 = xin1
    # du2 = uout1
    # dy1 = yin2
    # dv2 = vout2
    # du1 = uin3
    # dx2 = xout3
    # dv1 = vin4
    # dy2 = yout4

    # # Parameters with value
    # # r1) 
    # # dx1 = xin1
    # # du2 = uout1
    # # r2) 
    # # dy1 = yin2
    # # dv2 = vout2
    # # r3) 
    # # du1 = uin3
    # # dx2 = xout3
    # # r4) 
    # # dv1 = vin4
    # # dy2 = yout4

    # What rays to use where?
    # Ray 1 and 2 interact with A & C
    # Ray 3 and 4 interact with B & D

    print('xin = ',xin1)
    print('xout = ',xout1)
    print('yin = ',yin2)
    print('yout = ',yout2)
    print('uin = ',uin3)
    print('uout = ',uout3)
    print('vin = ',vin4)
    print('vout = ',vout4)

    # Construct the differential matrix
    dMat = np.array([[xout1/xin1,xout2/yin2,xout3/uin3,xout4/vin4],
                     [yout1/xin1,yout2/yin2,yout3/uin3,yout4/vin4],
                     [uout1/xin1,uout2/yin2,uout3/uin3,uout4/vin4],
                     [vout1/xin1,vout2/yin2,vout3/uin3,vout4/vin4]])

    return dMat

def _grab_ray_data(ray):
    return ray[0],ray[1],ray[2],ray[3]


def construct_rtm(rayin0,rayout0,
                  rayin1,rayout1,
                  rayin2,rayout2,
                  rayin3,rayout3,
                  rayin4,rayout4):

    rayin = np.array([rayin1-rayin0,
                      rayin2-rayin0,
                      rayin3-rayin0,
                      rayin4-rayin0])

    rayin = np.transpose(rayin)
    # rayin = np.transpose(rayin)
    rayinrank = np.linalg.matrix_rank(rayin)

    # print('Ray In Matrix',rayin)

    if rayinrank < 4:
        # print('Matrix rank is less than 4')
        return np.zeros([4,4])

    # print('Rayin shape = ',rayin.shape)
    # print(rayin)
    
    rayout = np.array([rayout1-rayout0,
                       rayout2-rayout0,
                       rayout3-rayout0,
                       rayout4-rayout0])

    # if np.linalg.matrix_rank(rayin) == 4:
    #     print('Ray In Matrix')
    #     print(rayin)
    # if np.linalg.matrix_rank(rayout) == 4:
    #     print('Ray Out Matrix')
    #     print(rayout)
    
    rtm = np.matmul(rayout,np.linalg.inv(rayin)) 
    # print(rayout)
    
    return rtm
    
def InspectRays(rayset):

    import matplotlib.pyplot as plt

    plt.figure(figsize=[10,5])
    plt.subplot(121)
    plt.title('Position')
    plt.scatter(rayset[0,:],rayset[1,:])
    plt.subplot(122)
    plt.title('Angle')
    plt.scatter(rayset[2,:],rayset[3,:])
    plt.show()
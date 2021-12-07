import matplotlib.pyplot as plt
import numpy as np

x = np.loadtxt('parabolatest.zmx_50_z9_xray.txt')
y = np.loadtxt('parabolatest.zmx_50_z9_yray.txt')
z = np.loadtxt('parabolatest.zmx_50_z9_zray.txt')

# a = np.loadtxt('parabolatest.zmx_50_aray.txt')
# b = np.loadtxt('parabolatest.zmx_50_bray.txt')

size = 1e-2

plt.figure(figsize=[7,7])
plt.subplot(121)
plt.title('Image Location [x,y]')
plt.xlim([-size,size])
plt.ylim([-size,size])
plt.scatter(x,y)
plt.subplot(122)
plt.title('Image Angle [x,y]')
plt.scatter(x/z,y/z)
plt.xlim([-size,size])
plt.ylim([-size,size])
plt.show()
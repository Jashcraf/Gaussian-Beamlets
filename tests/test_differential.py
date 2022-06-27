import numpy as np
import glets.differential as dif


pth = '/Users/jashcraft/Desktop/gbd-data/refract/refract_test_ray'

ray1infn = pth+'1_data_in.txt'
ray1outfn = pth+'1_data.txt'

ray2infn = pth+'2_data_in.txt'
ray2outfn = pth+'2_data.txt'

ray3infn = pth+'3_data_in.txt'
ray3outfn = pth+'3_data.txt'

ray4infn = pth+'4_data_in.txt'
ray4outfn = pth+'4_data.txt'

ray0infn = pth+'0_data_in.txt'
ray0outfn = pth+'0_data.txt'

ray1in = dif.readrays(ray1infn)
ray1out = dif.readrays(ray1outfn)

ray2in = dif.readrays(ray2infn)
ray2out = dif.readrays(ray2outfn)

ray3in = dif.readrays(ray3infn)
ray3out = dif.readrays(ray3outfn)

ray4in = dif.readrays(ray4infn)
ray4out = dif.readrays(ray4outfn)

ray0in = dif.readrays(ray0infn)
ray0out = dif.readrays(ray0outfn)

# dif.InspectRays(ray1in)
# dif.InspectRays(ray1out)
# dif.InspectRays(ray2in)
# dif.InspectRays(ray2out)
# dif.InspectRays(ray3in)
# dif.InspectRays(ray3out)
# dif.InspectRays(ray4in)
# dif.InspectRays(ray4out)

dif.InspectRays(ray1in)
dif.InspectRteays(ray1out)
dif.InspectRays(ray1in-ray0in)
dif.InspectRays(ray1out-ray0out)



# print(ray4in.shape)
raydata = ray4in
x = raydata[0,:]
y = raydata[1,:]
a = raydata[2,:]
b = raydata[3,:]
coords = ray4in[:,(x**2 + y**2 + a**2 + b**2)**(1/2) > 0] # ignores gut/invalid rays

RTM = np.zeros([4,4,ray4out.shape[1]])
import time
for i in range(int(ray4out.shape[1]/4)): # ,ray1in.shape[1]

    RTM[:,:,i] = dif.ConstructDiffMat(
                            ray0in[:,i],ray0out[:,i],
                            ray1in[:,i],ray1out[:,i],
                            ray2in[:,i],ray2out[:,i],
                            ray3in[:,i],ray3out[:,i],
                            ray4in[:,i],ray4out[:,i])

    # time.sleep(.5)
    print(RTM[:,:,i])

import glets.utils as ut

# print(RTM[:,:,400])

# ut.fourbyfour(RTM,coords=coords,size=6.5)
function rays = convert_ray_data(xData,yData,zData,lData,mData,nData)

x = xData;
y = yData;
z = zData; % unused
u = tan(pi/2 - acos(lData));
v = tan(pi/2 - acos(mData));
w = nData; % unused

rays = [x;y;u;v];

end
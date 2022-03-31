function ABCD = compute_differential(rayin0,rayout0,rayin1,rayout1,rayin2,rayout2,rayin3,rayout3,rayin4,rayout4)

disp('give me a breakpoint here')

% Compute Finite Difference
rayin1 = rayin1 - rayin0;
rayin2 = rayin2 - rayin0;
rayin3 = rayin3 - rayin0;
rayin4 = rayin4 - rayin0;

rayout1 = rayout1 - rayout0;
rayout2 = rayout2 - rayout0;
rayout3 = rayout3 - rayout0;
rayout4 = rayout4 - rayout0;

% Select ray data
[xin1,yin1,uin1,vin1] = grab_ray_data(rayin1);
[xin2,yin2,uin2,vin2] = grab_ray_data(rayin2);
[xin3,yin3,uin3,vin3] = grab_ray_data(rayin3);
[xin4,yin4,uin4,vin4] = grab_ray_data(rayin4);

[xout1,yout1,uout1,vout1] = grab_ray_data(rayout1);
[xout2,yout2,uout2,vout2] = grab_ray_data(rayout2);
[xout3,yout3,uout3,vout3] = grab_ray_data(rayout3);
[xout4,yout4,uout4,vout4] = grab_ray_data(rayout4);

% compute differential
ABCD = [xout1./xin1,xout2./yin2,xout3./uin3,xout4./vin4;...
        yout1./xin1,yout2./yin2,yout3./uin3,yout4./vin4;...
        uout1./xin1,uout2./yin2,uout3./uin3,uout4./vin4;...
        vout1./xin1,vout2./yin2,vout3./uin3,vout4./vin4];

end

function [xData,yData,uData,vData] = grab_ray_data(ray)

xData = ray(1,:);
yData = ray(2,:);
uData = ray(3,:);
vData = ray(4,:);

end
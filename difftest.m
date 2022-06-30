ray = [1,1,0,0];
dspace = logspace(1e-10,1,10);

Axx = [];
Axy = [];
Ayx = [];
Ayy = [];
Bxx = [];
Bxy = [];
Byx = [];
Byy = [];
Cxx = [];
Cxy = [];
Cyx = [];
Cyy = [];
Dxx = [];
Dxy = [];
Dyx = [];
Dyy = [];

for ijk = 1:length(dspace)
    
d = dspace(ijk);

wx = [1+d,1,0,0];
wy = [1,1+d,0,0];
dx = [1,1,d,0];
dy = [1,1,0,d];

efl = 25.4e-3;
dis  = efl;

sys = [0,0,dis,0;0,0,0,dis;-1/efl,0,1,0;0,-1/efl,0,1];

ray_prop = sys * ray';
wx_prop = sys * wx' - ray_prop;
wy_prop = sys * wy' - ray_prop;
dx_prop = sys * dx' - ray_prop;
dy_prop = sys * dy' - ray_prop;

xin1 = wx(1) - ray(1);
yin1 = wx(2) - ray(2);
uin1 = wx(3) - ray(3);
vin1 = wx(4) - ray(4);

xin2 = wy(1) - ray(1);
yin2 = wy(2) - ray(2);
uin2 = wy(3) - ray(3);
vin2 = wy(4) - ray(4);

xin3 = dx(1) - ray(1);
yin3 = dx(2) - ray(2);
uin3 = dx(3) - ray(3);
vin3 = dx(4) - ray(4);

xin4 = dy(1) - ray(1);
yin4 = dy(2) - ray(2);
uin4 = dy(3) - ray(3);
vin4 = dy(4) - ray(4);

xout1 = wx_prop(1);
yout1 = wx_prop(2);
uout1 = wx_prop(3);
vout1 = wx_prop(4);

xout2 = wy_prop(1);
yout2 = wy_prop(2);
uout2 = wy_prop(3);
vout2 = wy_prop(4);

xout3 = dx_prop(1);
yout3 = dx_prop(2);
uout3 = dx_prop(3);
vout3 = dx_prop(4);

xout4 = dy_prop(1);
yout4 = dy_prop(2);
uout4 = dy_prop(3);
vout4 = dy_prop(4);

abcd = [xout1/xin1,xout2/yin2,xout3/uin3,xout4/vin4;...
     yout1/xin1,yout2/yin2,yout3/uin3,yout4/vin4;...
     uout1/xin1,uout2/yin2,uout3/uin3,uout4/vin4;...
     vout1/xin1,vout2/yin2,vout3/uin3,vout4/vin4];

diff = abcd-sys;

Axx(end+1) = diff(1,1);
Axy(end+1) = diff(1,2);
Ayx(end+1) = diff(2,1);
Ayy(end+1) = diff(2,2);

Bxx(end+1) = diff(1,3);
Bxy(end+1) = diff(1,4);
Byx(end+1) = diff(2,3);
Byy(end+1) = diff(2,4);

Cxx(end+1) = diff(3,1);
Cxy(end+1) = diff(3,2);
Cyx(end+1) = diff(4,1);
Cyy(end+1) = diff(4,2);

Dxx(end+1) = diff(3,3);
Dxy(end+1) = diff(3,4);
Dyx(end+1) = diff(4,3);
Dyy(end+1) = diff(4,4);

end

figure(1)
hold on
semilogx(dspace,Bxx)
semilogx(dspace,Bxy)
semilogx(dspace,Byx)
semilogx(dspace,Byy)

semilogx(dspace,Cxx)
semilogx(dspace,Cxy)
semilogx(dspace,Cyx)
semilogx(dspace,Cyy)

semilogx(dspace,Dxx)
semilogx(dspace,Dxy)
semilogx(dspace,Dyx)
semilogx(dspace,Dyy)
hold off
ylim([-1e-14,1e-14])
legend('Bxx','Bxy','Byx','Byy','Cxx','Cxy','Cyx','Cyy','Dxx','Dxy','Dyx','Dyy')
set(gca, 'XScale', 'log')






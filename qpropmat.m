%% Qpropogation of Gaussian Beamlets

% Gauussian Parameters
wl = 1e-6;
npts = 4; % Used to autopopulate our pupil function
wo = 3*wl;
th = wl/(pi*wo);

% System Parameters
f = 0.5;
d = 0.5; % Distance between Source & Detector
dia = 1; % of pupil [m]

% Pupil Parameters
Rin = .005;
Rout = .01; 

%% Generate a Circular Pupil
x = linspace(-.01,.01,100); % 10mm
y = linspace(-.01,.01,100);

[x,y] = meshgrid(x,y);

pupil = zeros(size(x));
for yind = 1:length(pupil(:,1))
   for xind = 1:length(pupil(1,:)) 
        if (sqrt(x(yind,xind).^2 + y(yind,xind).^2) <= Rout) && (sqrt(x(yind,xind).^2 + y(yind,xind).^2)>= Rin)
            pupil(yind,xind) = 1;
        end
   end
end
surf(x,y,pupil,'linestyle','none')
view([0 90])

%% Define Sampling Area - Rectangular Distribution
u = linspace(-.01,.01,100);
v = linspace(-.01,.01,100);
[u,v] = meshgrid(u,v);

xbeams = length(x)/npts;
ybeams = length(y)/npts;
    



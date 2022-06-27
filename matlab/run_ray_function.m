%% Run a Matlab Raytrace - User inputs
pth = 'C:/Users/UASAL-OPTICS/Desktop/Gaussian-Beamlets/'; 
fn_ini = 'refract_test';
fov_to_trace_x = 0; % degrees
fov_to_trace_y = 0; % degrees
where = 'output';
% where = 'output';
maxsurface = 2;
fov_max = 1; % degrees
nrays = 21;
dz = 1e-4;

%% Some translation
fn = strcat(pth,fn_ini);
fnex = strcat(fn,'.zmx');
if strcmp(where,'input')==1
   surface = 1;
elseif strcmp(where,'output')==1
    surface = maxsurface;
end

%% Set up the Differential Vectors

% Geometrically orthogonal
dPx = [0,1,0,0,0]*dz/3250;
dPy = [0,0,1,0,0]*dz/3250;
dHx = [0,0,0,1,0]*dz;
dHy = [0,0,0,0,1]*dz;

% ijk = 1 should trace the base ray
% ijk = 2 should trace the 

for ijk = 1:5
    
    hx = (fov_to_trace_x + dHx(ijk))/fov_max;
    hy = (fov_to_trace_y + dHy(ijk))/fov_max;
    [xData,yData,zData,lData,mData,nData,l2Data,m2Data,n2Data,n1,n2 ] = MATLAB_BatchRayTrace_ReadRayData(...
        surface,nrays,hx,hy,dPx(ijk),dPy(ijk),fnex);
    
    xData = xData';
    yData = yData';
    zData = zData';

    lData = lData';
    mData = mData';
    nData = nData';

    l2Data = l2Data';
    m2Data = m2Data';
    n2Data = n2Data';
    % Create a table with the data and variable names
    T = table(xData,yData,zData,lData,mData,nData,l2Data,m2Data,n2Data);%,...
    %     'ColumnNames',...
    %     { 'xData','yData','zData','lData','mData','nData','l2Data','m2Data','n2Data','n1','n2'} );

    % Write data to text file
    fn = 'C:/Users/UASAL-OPTICS/Desktop/gbd-data/';
    fn = strcat(fn,fn_ini);
    if strcmp(where,'input')==1
        
        full_filename = strcat(fn,sprintf('_ray%d_data_in.txt',ijk-1));
        
    elseif strcmp(where,'output')==1
        
        full_filename = strcat(fn,sprintf('_ray%d_data.txt',ijk-1));
        
    end
    writetable(T, full_filename)
    
    
    
end

%% Write the data to a text file that python can interpret


%% Show the Intercept on surface

% figure(1)
% subplot(131)
%     hold on
%     title('Ray Intercepts')
%     scatter(xData,yData)
%     xlabel('distance [mm]')
%     hold off
% subplot(132)
%     hold on
%     title('Ray Angles')
%     scatter(lData,mData)
%     xlabel('Slope')
%     hold off
% subplot(133)
%     hold on
%     title('Ray Intercepts')
%     scatter(l2Data,m2Data)
%     xlabel('Slope')
%     hold off
    
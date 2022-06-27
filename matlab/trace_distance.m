function ABCD = trace_distance(surface,nrays,hx,hy,dPx,dPy,fnex)

% A cell to store ray data
% 1 = parabasal ray
% 2 = +X ray
% 3 = +Y ray
% 4 = +U ray
% 5 = +V ray

rayin = {};
rayout = {};
ABCD = {};

for abc = 1:5
    % Trace before the desired surface
    [xData1,yData1,zData1,lData1,mData1,nData1,~,~,~,~,~ ] = MATLAB_BatchRayTrace_ReadRayData(...
        surface-1,nrays,hx(abc),hy(abc),dPx(abc),dPy(abc),fnex);
    rays = convert_ray_data(xData1,yData1,zData1,lData1,mData1,nData1);
    rayin{end+1} = rays;
end

for def = 1:5
    % Trace to the desired surface
    [xData2,yData2,zData2,lData2,mData2,nData2,~,~,~,~,~ ] = MATLAB_BatchRayTrace_ReadRayData(...
            surface,nrays,hx(def),hy(def),dPx(def),dPy(def),fnex);
    rays = convert_ray_data(xData2,yData2,zData2,lData2,mData2,nData2);
    rayout{end+1} = rays;
end

% Create ABCD Matrix
ray0in = rayin{1};
ray0out = rayout{1};
ray1in = rayin{2};
ray1out = rayout{2};
ray2in = rayin{3};
ray2out = rayout{3};
ray3in = rayin{4};
ray3out = rayout{4};
ray4in = rayin{5};
ray4out = rayout{5};

% Compute differential
    for ghi = 1:length(ray0in)
        ABCD{end+1} = compute_differential_distance(ray0in(:,ghi),ray0out(:,ghi),...
                                                    ray1in(:,ghi),ray1out(:,ghi),...
                                                    ray2in(:,ghi),ray2out(:,ghi),...
                                                    ray3in(:,ghi),ray3out(:,ghi),...
                                                    ray4in(:,ghi),ray4out(:,ghi));
          
    end

end
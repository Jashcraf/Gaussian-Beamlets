function ABCD = trace_refraction(surface,nrays,hx,hy,dPx,dPy,fnex)

% A cell to store ray data
% 1 = parabasal ray
% 2 = +X ray
% 3 = +Y ray
% 4 = +U ray
% 5 = +V ray
ABCD = {};
rayin = {};
rayout = {};

    for abc = 1:5
       [xData,yData,zData,lData,mData,nData,l2Data,m2Data,n2Data,~,~ ] = MATLAB_BatchRayTrace_ReadRayData(...
            surface-1,nrays,hx(abc),hy(abc),dPx(abc),dPy(abc),fnex); 

        % Compute AOI from AOE
        numerator = (lData.*l2Data + mData.*m2Data + nData.*n2Data);
        denominator = sqrt(lData.^2 + mData.^2 + nData.^2).*sqrt(l2Data.^2 + m2Data.^2 + n2Data.^2);
        aoe = numerator./denominator;
        aoe = aoe - (aoe(1:length(xData)) > pi/2) * pi;
        aoe = abs(aoe);
        aoi = -aoe; % Snells law for reflection

        % Derive K Vectors
        kout = [lData;mData;nData]; % ray out
        norm = -[l2Data;m2Data;n2Data];
        kin = kout - 2*cos(aoi).*norm;

        % Convert to rays
        rayin{end+1} = convert_ray_data(xData,yData,zData,kin(1,:),kin(2,:),kin(3,:));
        rayout{end+1} = convert_ray_data(xData,yData,zData,lData,mData,nData);


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

% Compute differential in a for loop
% length() returns largest dimension of the array you pass it
    for def = 1:length(ray0in)
        ABCD{end+1} = compute_differential_refraction(ray0in(:,def),ray0out(:,def),...
                                                      ray1in(:,def),ray1out(:,def),...
                                                      ray2in(:,def),ray2out(:,def),...
                                                      ray3in(:,def),ray3out(:,def),...
                                                      ray4in(:,def),ray4out(:,def));
    end

end
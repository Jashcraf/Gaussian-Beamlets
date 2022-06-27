function [ xData,yData,zData,lData,mData,nData,l2Data,m2Data,n2Data,n1,n2,Px,Py ] = MATLAB_BatchRayTrace_ReadRayData( surface, nrays, hx, hy, dPx, dPy, filename )
% This is just to initialize the ZOS-API connection
if ~exist('args', 'var')
    args = [];
end

% Initialize the OpticStudio connection
TheApplication = InitConnection();
if isempty(TheApplication)
    % failed to initialize a connection
    r = [];
else
    try
        [ xData,yData,zData,lData,mData,nData,l2Data,m2Data,n2Data,n1,n2,Px,Py ] = BeginApplication(TheApplication, surface, nrays, hx, hy, dPx, dPy, filename);
        CleanupConnection(TheApplication);
    catch err
        CleanupConnection(TheApplication);
        rethrow(err);
    end
end
end




function [ xData,yData,zData,lData,mData,nData,l2Data,m2Data,n2Data,n1,n2,Px,Py ] = BeginApplication(TheApplication, surface, nrays, hx, hy, dPx, dPy, filename)

    import ZOSAPI.*;
    
    % Creates empty system
    TheSystem = TheApplication.CreateNewSystem(ZOSAPI.SystemType.Sequential);

    % Add your custom code here...
    
    % opens desired ZMX file
    % filename = '/Double Gauss 28 degree field.zmx';
    % TheSystem.LoadFile(System.String.Concat(TheApplication.SamplesDir, filename), false);
    TheSystem.LoadFile(filename,false)
    disp('File to load')
    disp(filename)
    
    
    % Set up Batch Ray Trace
    % creates batch raytrace in API
	raytrace = TheSystem.Tools.OpenBatchRayTrace();
    nsur = TheSystem.LDE.NumberOfSurfaces; % Ignore the object surface?
    if surface >= nsur
        disp('Final surface exceeds number of surfaces')
        disp('Bypassing by tracing to surface:')
        surface = nsur;
        disp(surface)
    end
    
    max_rays = nrays;
    total_rays_in_both_axes = (max_rays) * (max_rays);
    
    % creates batch raytrace in API
    % Performs a batch unpolarized ray trace, using normalized pupil coordiantes; this is similar to the DDE ray trace command, mode 0.
    % CreateNormUnpol (int MaxRays, RaysType rayType, int toSurface) what
    % are rayTypes?
    RayTraceData = raytrace.CreateNormUnpol(total_rays_in_both_axes, ZOSAPI.Tools.RayTrace.RaysType.Real, surface);
    
	% offloads processing to C# dll 
    NET.addAssembly(strcat(pwd, '\RayTrace.dll')); % Where does NET come from?
    import BatchRayTrace.*;
    
    % Attach the helper class to the ray database reader
    tic(); % Begin timing
    close all;
    % color_ary = {'blue', 'green', 'red', 'gold', 'pink', 'cyan', 'purple', 'teal'};
    
    
    wave = 1;
    % Call the function ReadNormUnpolData from dll
    % The syntax is ReadNormUnpolData(IBatchRayTrace rt,
    % IRayTraceNormUnpolData rtt) Unpol
    dataReader = ReadNormUnpolData(raytrace, RayTraceData); % Notice how it isn't in a for loop
    
    % The variable dataReader can now use all the functions defined in the DLL
    dataReader.ClearData();
    
    % On-axis field
    x_field = hx;
    y_field = hy;
    
    disp('Tracing field')
    
    disp(x_field)
    disp(y_field)

    % creates array of field (Hx/Hy) and pupil (Px/Py) coordinates
    Hx = ones(total_rays_in_both_axes, 1).*x_field;
    Hy = ones(total_rays_in_both_axes, 1).*y_field;
    
    % pupil fill factor
    pfac = 0.99;
    disp('pupil fill factor')
    disp(pfac)

    % square grid pattern
    Pxi = reshape(bsxfun(@times, linspace(1, 1, max_rays)', linspace(-1, 1, max_rays)), [max_rays^2, 1]);
    Pyi = reshape(bsxfun(@times, linspace(1, -1, max_rays)', linspace(1, 1, max_rays)), [max_rays^2, 1]);
    Pxi = Pxi*pfac;
    Pyi = Pyi*pfac;
    
    % Apply differentials
    Pxi = Pxi + dPx;
    Pyi = Pyi + dPy;
    
    % limits pupil to unit circle
    Px = Pxi(1:total_rays_in_both_axes) .* (sqrt((Pxi(1:total_rays_in_both_axes).^2) + Pyi(1:total_rays_in_both_axes).^2) <= 1);
    Py = Pyi(1:total_rays_in_both_axes) .* (sqrt((Pxi(1:total_rays_in_both_axes).^2) + Pyi(1:total_rays_in_both_axes).^2) <= 1);
    
    % for dithered grid, uncomment next line
    %[Px, Py] = rand_circle(total_rays_in_both_axes);

    % converts from matlab arrays to .NET arrays 
    HxNet = NET.convertArray(Hx, 'System.Double');
    HyNet = NET.convertArray(Hy, 'System.Double');
    PxNet = NET.convertArray(Px, 'System.Double');
    PyNet = NET.convertArray(Py, 'System.Double');

    dataReader.AddRay(wave, HxNet, HyNet, PxNet, PyNet, ZOSAPI.Tools.RayTrace.OPDMode.None);

    % Configure the maximum number of segments to read at one time.
    % Note that this is a tradeoff between speed and memory usage
    rayData = dataReader.InitializeOutput(max_rays);
    isFinished = false;
    totalRaysRead = 0;
    maxRays = 0;

    while ~isFinished
        % fill the next block of data
        readSegments = dataReader.ReadNextBlock(rayData);

        if readSegments == 0
            isFinished = true;
        else
            % Note - MATLAB arrays are 1-based, however .NET arrays are 0-based, so we have to be carefull...
            maxRays = max(rayData.rayNumber.double);
            
            % Position Data
            xData = rayData.X.double;
            yData = rayData.Y.double;
            zData = rayData.Z.double;
            
            % Grab the surfaces Rotation Matrix and Offset
%             [success,R11,R12,R13,R21,R22,R23,R31,R32,R33,XO,YO,ZO] = TheSystem.LDE.GetGlobalMatrix(surface);
%             if success ~= 1
%                disp('Global Matrix Failure') 
%             end
%             Rmat = [R11,R12,R13;R21,R22,R23;R31,R32,R33];
%             Dmat = [xData_in;yData_in;zData_in];
%             Offset = zeros(size(Dmat));
%             Offset(1,:) = XO;
%             Offset(2,:) = YO;
%             Offset(3,:) = ZO;
%             dOut = Offset + Rmat*Dmat; 
%             
%             xData = dOut(1,:);
%             yData = dOut(2,:);
%             zData = dOut(3,:);
            
            
            % Direction Cosine Data after requested surface
            lData = rayData.L.double;
            mData = rayData.M.double;
            nData = rayData.N.double;
            
            % vector normal of the ray at the requested surface. 
            l2Data = rayData.l2.double;
            m2Data = rayData.m2.double;
            n2Data = rayData.n2.double;
            
            % calculate index at surface
            % refractive indices
            n1 = TheSystem.MFE.GetOperandValue(ZOSAPI.Editors.MFE.MeritOperandType.INDX, surface - 1, wave, 0, 0, 0, 0, 0, 0);
            n2 = TheSystem.MFE.GetOperandValue(ZOSAPI.Editors.MFE.MeritOperandType.INDX, surface, wave, 0, 0, 0, 0, 0, 0);
            

        end

        if totalRaysRead >= maxRays
            isFinished = true;
        end    
    end
    
raytrace.Close();

toc();
disp(['Total Rays: ', num2str(maxRays)]);
%r = [xData,yData,zData,lData,mData,nData,l2Data,m2Data,n2Data];

end

function app = InitConnection()

import System.Reflection.*;

% Find the installed version of OpticStudio.
zemaxData = winqueryreg('HKEY_CURRENT_USER', 'Software\Zemax', 'ZemaxRoot');
NetHelper = strcat(zemaxData, '\ZOS-API\Libraries\ZOSAPI_NetHelper.dll');
% Note -- uncomment the following line to use a custom NetHelper path
% NetHelper = '@{NETHELP}';
NET.addAssembly(NetHelper);

success = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize();
% Note -- uncomment the following line to use a custom initialization path
% success = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize('C:\Program Files\Zemax OpticStudio\');
if success == 1
    LogMessage(strcat('Found OpticStudio at: ', char(ZOSAPI_NetHelper.ZOSAPI_Initializer.GetZemaxDirectory())));
else
    app = [];
    return;
end

% Now load the ZOS-API assemblies
NET.addAssembly(AssemblyName('ZOSAPI_Interfaces'));
NET.addAssembly(AssemblyName('ZOSAPI'));

% Create the initial connection class
TheConnection = ZOSAPI.ZOSAPI_Connection();

% Attempt to create a Standalone connection

% NOTE - if this fails with a message like 'Unable to load one or more of
% the requested types', it is usually caused by try to connect to a 32-bit
% version of OpticStudio from a 64-bit version of MATLAB (or vice-versa).
% This is an issue with how MATLAB interfaces with .NET, and the only
% current workaround is to use 32- or 64-bit versions of both applications.
app = TheConnection.CreateNewApplication();
if isempty(app)
   HandleError('An unknown connection error occurred!');
end
if ~app.IsValidLicenseForAPI
    HandleError('License check failed!');
    app = [];
end

end

function LogMessage(msg)
disp(msg);
end

function HandleError(error)
ME = MXException(error);
throw(ME);
end

function  CleanupConnection(TheApplication)
% Note - this will close down the connection.

% If you want to keep the application open, you should skip this step
% and store the instance somewhere instead.
TheApplication.CloseApplication();
end

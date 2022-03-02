import clr, os, winreg
from itertools import islice

import matplotlib.pyplot as plt
import numpy as np
from System import Enum, Int32, Double

class PythonStandaloneApplication(object):
    class LicenseException(Exception):
        pass
    class ConnectionException(Exception):
        pass
    class InitializationException(Exception):
        pass
    class SystemNotPresentException(Exception):
        pass

    def __init__(self, path=None):
        # determine location of ZOSAPI_NetHelper.dll & add as reference
        aKey = winreg.OpenKey(winreg.ConnectRegistry(None, winreg.HKEY_CURRENT_USER), r"Software\Zemax", 0, winreg.KEY_READ)
        zemaxData = winreg.QueryValueEx(aKey, 'ZemaxRoot')
        NetHelper = os.path.join(os.sep, zemaxData[0], r'ZOS-API\Libraries\ZOSAPI_NetHelper.dll')
        winreg.CloseKey(aKey)
        clr.AddReference(NetHelper)
        import ZOSAPI_NetHelper

        # Find the installed version of OpticStudio
        if path is None:
            isInitialized = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize()
        else:
            # Note -- uncomment the following line to use a custom initialization path
            isInitialized = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize(path)

        # determine the ZOS root directory
        if isInitialized:
            dir = ZOSAPI_NetHelper.ZOSAPI_Initializer.GetZemaxDirectory()
        else:
            raise PythonStandaloneApplication.InitializationException("Unable to locate Zemax OpticStudio.  Try using a hard-coded path.")

        # add ZOS-API referencecs
        clr.AddReference(os.path.join(os.sep, dir, "ZOSAPI.dll"))
        clr.AddReference(os.path.join(os.sep, dir, "ZOSAPI_Interfaces.dll"))
        import ZOSAPI

        # create a reference to the API namespace
        self.ZOSAPI = ZOSAPI

        # create a reference to the API namespace
        self.ZOSAPI = ZOSAPI

        # Create the initial connection class
        self.TheConnection = ZOSAPI.ZOSAPI_Connection()

        if self.TheConnection is None:
            raise PythonStandaloneApplication.ConnectionException("Unable to initialize .NET connection to ZOSAPI")

        self.TheApplication = self.TheConnection.CreateNewApplication()
        if self.TheApplication is None:
            raise PythonStandaloneApplication.InitializationException("Unable to acquire ZOSAPI application")

        if self.TheApplication.IsValidLicenseForAPI == False:
            raise PythonStandaloneApplication.LicenseException("License is not valid for ZOSAPI use")

        self.TheSystem = self.TheApplication.PrimarySystem
        if self.TheSystem is None:
            raise PythonStandaloneApplication.SystemNotPresentException("Unable to acquire Primary system")

    def __del__(self):
        if self.TheApplication is not None:
            self.TheApplication.CloseApplication()
            self.TheApplication = None

        self.TheConnection = None

    def OpenFile(self, filepath, saveIfNeeded):
        if self.TheSystem is None:
            raise PythonStandaloneApplication.SystemNotPresentException("Unable to acquire Primary system")
        self.TheSystem.LoadFile(filepath, saveIfNeeded)

    def CloseFile(self, save):
        if self.TheSystem is None:
            raise PythonStandaloneApplication.SystemNotPresentException("Unable to acquire Primary system")
        self.TheSystem.Close(save)

    def SamplesDir(self):
        if self.TheApplication is None:
            raise PythonStandaloneApplication.InitializationException("Unable to acquire ZOSAPI application")

        return self.TheApplication.SamplesDir

    def ExampleConstants(self):
        if self.TheApplication.LicenseStatus == self.ZOSAPI.LicenseStatusType.PremiumEdition:
            return "Premium"
        elif self.TheApplication.LicenseStatus == self.ZOSAPI.LicenseStatusTypeProfessionalEdition:
            return "Professional"
        elif self.TheApplication.LicenseStatus == self.ZOSAPI.LicenseStatusTypeStandardEdition:
            return "Standard"
        else:
            return "Invalid"

    def reshape(self, data, x, y, transpose = False):
        """Converts a System.Double[,] to a 2D list for plotting or post processing

        Parameters
        ----------
        data      : System.Double[,] data directly from ZOS-API
        x         : x width of new 2D list [use var.GetLength(0) for dimension]
        y         : y width of new 2D list [use var.GetLength(1) for dimension]
        transpose : transposes data; needed for some multi-dimensional line series data

        Returns
        -------
        res       : 2D list; can be directly used with Matplotlib or converted to
                    a numpy array using numpy.asarray(res)
        """
        if type(data) is not list:
            data = list(data)
        var_lst = [y] * x;
        it = iter(data)
        res = [list(islice(it, i)) for i in var_lst]
        if transpose:
            return self.transpose(res);
        return res

    def transpose(self, data):
        """Transposes a 2D list (Python3.x or greater).

        Useful for converting mutli-dimensional line series (i.e. FFT PSF)

        Parameters
        ----------
        data      : Python native list (if using System.Data[,] object reshape first)

        Returns
        -------
        res       : transposed 2D list
        """
        if type(data) is not list:
            data = list(data)
        return list(map(list, zip(*data)))


"""
Paste code there
"""

if __name__ == '__main__':
    zos = PythonStandaloneApplication()

    # User inputs to the beamlet propagation
    max_rays = 50
    filename = "parabolatest.zmx"

    # Just a single fov at a time
    fovx_deg = 0.0
    fovy_deg = 0.0

    # load local variables
    ZOSAPI = zos.ZOSAPI
    TheApplication = zos.TheApplication
    TheSystem = zos.TheSystem
    if not os.path.exists(TheApplication.SamplesDir + "\\API\\Python"):
        os.makedirs(TheApplication.SamplesDir + "\\API\\Python")

    # Set up primary optical system
    """
    WARNING THIS ONLY APPEARS TO WORK WHEN THE LENS FILE IS IN THE SAMPLES DIRECTORY UNDER SEQUENTIAL/OBJECTIVES/FILE.ZMX
    Why does this happen
    whyyyyyy
    """
    sampleDir = TheApplication.SamplesDir
    file = "parabolatest.zmx"
    testFile = sampleDir + "\\Sequential\\Objectives\\" + file
    TheSystem.LoadFile(testFile, False)
    # sampleDir = TheApplication.SamplesDir
    # file = filename
    # testFile = sampleDir + filename #sampleDir + "\\Sequential\\Objectives\\" + file
    # TheSystem.LoadFile(filename, False)

    #! [e22s01_py]
    # Set up Batch Ray Trace
    raytrace = TheSystem.Tools.OpenBatchRayTrace()
    nsur = TheSystem.LDE.NumberOfSurfaces
    print(nsur)
    normUnPolData = raytrace.CreateNormUnpol((max_rays + 1) * (max_rays + 1), ZOSAPI.Tools.RayTrace.RaysType.Real, nsur)
    #! [e22s01_py]

    #! [e22s02_py]
    # Define batch ray trace constants
    hx = fovx_deg
    max_wave = 1# TheSystem.SystemData.Wavelengths.NumberOfWavelengths
    num_fields = 1 #TheSystem.SystemData.Fields.NumberOfFields
    hy_ary = fovy_deg #np.array([0, 0.707, 1])
    #! [e22s02_py]
    waveNumber = 1
    """
    Don't think this is neccessary because we just want ray data
    """
    # Initialize x/y/z image plane arrays -
    num_fields = 1
    #max_wave = 0
    x_ary = np.empty((max_rays + 1) * (max_rays + 1))
    y_ary = np.empty((max_rays + 1) * (max_rays + 1))
    z_ary = np.empty((max_rays + 1) * (max_rays + 1))

    # Initialize alpha/beta/gamma image plane arrays -
    a_ary = np.empty((max_rays + 1) * (max_rays + 1))
    b_ary = np.empty((max_rays + 1) * (max_rays + 1))
    g_ary = np.empty((max_rays + 1) * (max_rays + 1))

    #! [e22s03_py]
    # Determine maximum field in Y only
    max_field = 0.0
    # for i in range(1, num_fields + 1):
    #     if (TheSystem.SystemData.Fields.GetField(i).Y > max_field):
    #         max_field = TheSystem.SystemData.Fields.GetField(i).Y
    #! [e22s03_py]

    # Determine number of rays in the entrance pupil - assume a square
    px = np.linspace(-1,1,(max_rays + 1))
    px,py = np.meshgrid(px,px)
    px = np.ravel(px)
    py = np.ravel(py)

    # plt.rcParams["figure.figsize"] = (15, 4)
    # colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')

    # Add rays in a for-loop based on pupil & image location
    normUnPolData.AddRay(waveNumber, hx, hy_ary, px[0], py[0], Enum.Parse(ZOSAPI.Tools.RayTrace.OPDMode, "None"))

    if TheSystem.SystemData.Fields.GetFieldType() == ZOSAPI.SystemData.FieldType.Angle:
        field_type = 'Angle'
    elif TheSystem.SystemData.Fields.GetFieldType() == ZOSAPI.SystemData.FieldType.ObjectHeight:
        field_type = 'Height'
    elif TheSystem.SystemData.Fields.GetFieldType() == ZOSAPI.SystemData.FieldType.ParaxialImageHeight:
        field_type = 'Height'
    elif TheSystem.SystemData.Fields.GetFieldType() == ZOSAPI.SystemData.FieldType.RealImageHeight:
        field_type = 'Height'

    #! [e22s04_py]
    # Adding Rays to Batch, varying normalised object height hy
    normUnPolData.ClearData()
    waveNumber = 0
    #for i = 1:((max_rays + 1) * (max_rays + 1))
    for i in range(1, (max_rays) * (max_rays) + 1):

        # px = np.random.random() * 2 - 1
        # py = np.random.random() * 2 - 1

        # while (px*px + py*py > 1):
        #     py = np.random.random() * 2 - 1
        normUnPolData.AddRay(waveNumber, hx, hy_ary, px[i], py[i], Enum.Parse(ZOSAPI.Tools.RayTrace.OPDMode, "None"))
    #! [e22s04_py]

    raytrace.RunAndWaitForCompletion()

    #! [e22s05_py]
    # Read batch raytrace and display results
    normUnPolData.StartReadingResults()

    # Python NET requires all arguments to be passed in as reference, so need to have placeholders
    wave = max_wave
    field = max_field
    sysInt = Int32(1)
    sysDbl = Double(1.0)
    # output structure is
    # 0) success
    # 1) ray number
    # 2) errCode
    # 3) vignetteCode
    # 4) x pos
    # 5) x ang
    # 6) y pos
    # 7) y ang
    # 8) z pos
    # 9) z ang
    # 10) intensity -> use this to get amps?
    output = normUnPolData.ReadNextResult(sysInt, sysInt, sysInt,
                   sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl)
    

    xlist = []
    ylist = []
    zlist = []

    while output[0]:                                                    # success
        if ((output[2] == 0) and (output[3] == 0)):                     # ErrorCode & vignetteCode
            # store in a given raynumber
            x_ary[output[1] - 1] = output[4]   # X
            y_ary[output[1] - 1] = output[5]   # Y
            z_ary[output[1] - 1] = output[6]   # Z
            a_ary[output[1] - 1] = output[7]   # Alpha
            b_ary[output[1] - 1] = output[8]   # Beta
            g_ary[output[1] - 1] = output[9]   # Gamma
            xlist.append(output[4])
            ylist.append(output[5])
            zlist.append(output[6])

        output = normUnPolData.ReadNextResult(sysInt, sysInt, sysInt,
                   sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl)
            #! [e22s05_py]
            # temp = plt.plot(np.squeeze(x_ary[field - 1, wave - 1, :]), np.squeeze(y_ary[field - 1, wave - 1, :]), '.', ms = 1, color = colors[wave - 1])
    
    np.savetxt('{fname}_{nrays}_z9_xray.txt'.format(fname=filename,nrays=max_rays),x_ary)
    np.savetxt('{fname}_{nrays}_z9_yray.txt'.format(fname=filename,nrays=max_rays),y_ary)
    np.savetxt('{fname}_{nrays}_z9_zray.txt'.format(fname=filename,nrays=max_rays),z_ary)
    np.savetxt('{fname}_{nrays}_z9_aray.txt'.format(fname=filename,nrays=max_rays),a_ary)
    np.savetxt('{fname}_{nrays}_z9_bray.txt'.format(fname=filename,nrays=max_rays),b_ary)
    np.savetxt('{fname}_{nrays}_z9_gray.txt'.format(fname=filename,nrays=max_rays),g_ary)
    del zos
    zos = None
    # This isn't in a function rn
    # return x_ary,y_ary,z_ary,a_ary,b_ary,g_ary

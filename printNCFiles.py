from netCDF4 import Dataset
import numpy as np
def createNC(outFile,T,Y,X):
    NY = np.size(Y)
    NX = np.size(X)
    NT = np.size(T)
    newf = Dataset(outFile, 'w', format='NETCDF3_64BIT')
# CREATE DIMENSIONS
    time = newf.createDimension('time',None)
    y = newf.createDimension('y',NY)
    x = newf.createDimension('x',NX)
# CREATE DIMENSION VARIABLES
    y = newf.createVariable('y','f4',('y'))
    x = newf.createVariable('x','f4',('x'))
    times = newf.createVariable('time','f4',('time'))
# CREATE VARIABLES
    v = newf.createVariable('v','f4',('time','y','x'))
    u = newf.createVariable('u','f4',('time','y','x'))
    vvar = newf.createVariable('vvar','f4',('time','y','x'))
    uvar = newf.createVariable('uvar','f4',('time','y','x'))
# INITIALIZE VARIABLES
    y[:] = np.squeeze(Y)
    x[:] = np.squeeze(X)
    times[:] = np.squeeze(T) 
    newf.close()

# write on netcdf file ###########################################################
def writeNC(f,varname,data):
   NT = np.size(data,0)
   var = f.variables[varname]
   var[0:NT,:,:,:] = data
   return f

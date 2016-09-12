import sys
sys.path.append('/nethome/swang/python/lib/python2.7/site-packages')
from netCDF4 import Dataset
import numpy as np
def createNC(outFile,T,Y,X,hyp_v,hyp_u,hyp_names):
    NY = np.size(Y)
    NX = np.size(X)
    NT = np.size(T)
    NH = np.size(hyp_v)
    NHn= np.size(hyp_names)
    newf = Dataset(outFile, 'w', format='NETCDF3_64BIT')
# CREATE DIMENSIONS
    time = newf.createDimension('time',None)
    y = newf.createDimension('y',NY)
    x = newf.createDimension('x',NX)
    h = newf.createDimension('hyperparam',NH)
    hn = newf.createDimension('param_names',NHn)
# CREATE DIMENSION VARIABLES
    y = newf.createVariable('y','f4',('y'))
    x = newf.createVariable('x','f4',('x'))
    times = newf.createVariable('time','f4',('time'))
    hu = newf.createVariable('hyperparam_u','f4',('hyperparam'))
    hv = newf.createVariable('hyperparam_v','f4',('hyperparam'))
    hn = newf.createVariable('param_names','f4',('param_names'))
# CREATE VARIABLES
    v = newf.createVariable('v','f4',('time','y','x'))
    u = newf.createVariable('u','f4',('time','y','x'))
    vvar = newf.createVariable('vvar','f4',('time','y','x'))
    uvar = newf.createVariable('uvar','f4',('time','y','x'))
# INITIALIZE VARIABLES
    y[:] = np.squeeze(Y)
    x[:] = np.squeeze(X)
    times[:] = np.squeeze(T) 
    hv[:] = np.squeeze(hyp_v)
    hu[:] = np.squeeze(hyp_u)
    hn[:] = np.squeeze(hyp_names)
    newf.close()

# write on netcdf file ###########################################################
def writeNC(f,varname,data):
   NT = np.size(data,0)
   var = f.variables[varname]
   var[0:NT,:,:] = data
   return f

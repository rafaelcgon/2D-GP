import numpy as np
from datetime import datetime
from laser_class import drifter,interpolated_tracks
from laser_io_methods import save_object
import cPickle as pickle
import os
from dateutil.parser import parse
from datetime import datetime,timedelta
from subprocess import Popen, PIPE
from scipy import interpolate
import GPy
# open files
filename = 'ALL_2016_2_7.pkl'
kriging_filename = 'kriging2_' + filename
spline_filename = 'spline_' + filename
with open(kriging_filename,'rb') as input:
            kr = pickle.load(input)

with open(spline_filename,'rb') as input:
            spl = pickle.load(input)

# compare interpolation methods
err_lon = np.square(kr.lon-spl.lon)
err_lat = np.square(kr.lat-spl.lat)
err_lon[np.where(kr.drogueStat==0)]=np.nan
err_lat[np.where(kr.drogueStat==0)]=np.nan
RMSE_Lon = np.sqrt(np.nanmean(err_lon,axis=0))
RMSE_Lat = np.sqrt(np.nanmean(err_lat,axis=0))

# compare velocities
spl.u[np.where(kr.drogueStat==0)]=np.nan
spl.v[np.where(kr.drogueStat==0)]=np.nan
spl_u = (spl.u[:,1:]+spl.u[:,:-1])/2.
spl_v = (spl.v[:,1:]+spl.v[:,:-1])/2.
err_u = np.square(kr.u - spl_u)
err_v = np.square(kr.v-spl_v)
#err_u[np.where(kr.drogueStat==0)]=np.nan
#err_v[np.where(kr.drogueStat==0)]=np.nan
RMSE_u = np.sqrt(np.nanmean(err_u,axis=0))
RMSE_v = np.sqrt(np.nanmean(err_v,axis=0))

rmse_lon=np.zeros(np.size(kr.id))
rmse_lat=np.zeros(np.size(kr.id))
for i in range(np.size(kr.id)):   
    error_lo = np.square(kr.lon[0,:]-spl.lon[0,:])
    error_la = np.square(kr.lat[0,:]-spl.lat[0,:])
    rmse_lon[i] = np.sqrt(np.nanmean(error_lo))
    rmse_lat[i] = np.sqrt(np.nanmean(error_la))

# filter bad drifters out (mark as undrogued for now...)

ind = np.where(RMSE_Lon>0.001)[0]
for i in ind:
   bad = np.where(err_lon[:,i]>0.0001)[0]
   for j in bad:
# check if error persists in the next time step and if drifter is drogued
      if i < np.size(err_lon,1)-1:
         if (err_lon[j,i+1]>0.0001)&(kr.drogueStat[j,i] == 1): 
            print kr.drogueStat[j,i]
            kr.drogueStat[j,i:] = 0
            spl.drogueStat[j,i:] = 0
            print i,j,kr.id[j],kr.drogueStat[j,i]

save_object(kr,'new_'+kriging_filename)
save_object(spl,'new_'+spline_filename)



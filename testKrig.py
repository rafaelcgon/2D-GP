import numpy as np
from datetime import datetime
import cPickle as pickle
import GPy
import pyproj
from matplotlib import rc,rcParams
import myKernel
import scipy.io as sio
from printNCFiles import createNC,writeNC
import sys
from netCDF4 import Dataset
import memory_profiler as mprof
from krig import getData
import gc
from sklearn.gaussian_process import kernels,GaussianProcessRegressor


def kriging(st=0,et=48,sample_step=-2,skip=3,num_restarts=20,save=1,output='testModel',gpy=0):
   """
   Estimate velocity field using rbf kernels on t,y,x (K = Kt*Ky*Kx)
   st,et = initial,final time step of drifter data
   sample_step = pick data each sample_step. If negative, pick data at 
each time step = sample_step*(-1)     
   skip = skip drifters. If skip ==1, use all drifters
   num_restart: number of restarts on optimization
   save: if == 0, return values, if not, save as output+'.mat' and output+'.pkl' 
   """
#   st = 0; et = 144
   startTime = datetime.now()
   time,latt,lont,vob,uob,validPoints = getData(st,et)

   # origin of cartesian coord.
   lat0 = 28.8
   lon0 = -88.55
   NAD83=pyproj.Proj("+init=EPSG:3452") #Louisiana South (ftUS)
   xob,yob=NAD83(lont,latt)
   tob = time[:,None]; tob = np.repeat(tob,np.size(latt,1),axis=1)
   yob[np.where(np.isnan(lont))]=np.nan
   xob[np.where(np.isnan(lont))]=np.nan

   x_ori = np.nanmin(xob)+2; y_ori = np.nanmin(yob)+2
   xob = (xob - x_ori)/1000. # in km
   yob = (yob - y_ori)/1000. # in km

   Nd = np.size(xob,1)/skip # number of drifters
   if (sample_step < 0)|(skip>1): # get samples per time step
      ss = np.abs(sample_step)
      sample_step = 1
      samples = np.arange(0,np.size(tob,0),ss)
 
      if (skip>1): # !!!!
         testt = np.arange(np.size(tob,0)) # 
         testd = set(np.arange(0,np.size(tob,1))) - set(np.arange(0,np.size(tob,1),skip))
         testd = np.array(list(testd))
      else:
         testd = np.arange(0,np.size(tob,1)) 
         if ss > 1:
            testt = set(np.arange(0,np.size(tob,0))) - set(samples)
            testt = np.array(list(testt))
         else:
            testt = samples
      to = np.reshape(tob[samples,::skip],[-1,1])
      yo = np.reshape(yob[samples,::skip],[-1,1])
      xo = np.reshape(xob[samples,::skip],[-1,1])
      uo = np.reshape(uob[samples,::skip],[-1,1])
      vo = np.reshape(vob[samples,::skip],[-1,1])
      tt = np.reshape(tob[testt[:,None],testd],[-1,1])
      yt = np.reshape(yob[testt[:,None],testd],[-1,1])
      xt = np.reshape(xob[testt[:,None],testd],[-1,1])
      ut = np.reshape(uob[testt[:,None],testd],[-1,1])
      vt = np.reshape(vob[testt[:,None],testd],[-1,1])

   else:
      ss = 0
      samples = np.arange(0,xob.size,sample_step) #np.random.randint(0,xo.size,nsamples)
      test = set(np.arange(xob.size)) - set(samples)
      test = np.array(list(test))
      xt = np.reshape(xob,[-1])[test,None]
      yt = np.reshape(yob,[-1])[test,None]
      tt = np.reshape(tob,[-1])[test,None]
      ut = np.reshape(uob,[-1])[test,None]
      vt = np.reshape(vob,[-1])[test,None]
      xo = np.reshape(xob,[-1])[samples,None]
      yo = np.reshape(yob,[-1])[samples,None]
      to = np.reshape(tob,[-1])[samples,None]
      uo = np.reshape(uob,[-1])[samples,None] 
      vo = np.reshape(vob,[-1])[samples,None]

  
   validPoints = np.where((~np.isnan(xo))&(~np.isnan(yo)))
   to = to[validPoints][:,None]
   xo = xo[validPoints][:,None]
   yo = yo[validPoints][:,None]
   uo = uo[validPoints][:,None]
   vo = vo[validPoints][:,None]
   validPoints = np.where((~np.isnan(xt))&(~np.isnan(yt)))
   tt = tt[validPoints][:,None]
   xt = xt[validPoints][:,None]
   yt = yt[validPoints][:,None]
   ut = ut[validPoints][:,None]
   vt = vt[validPoints][:,None]

   print 'number of observations: '+str(np.size(vo))
   output_mat = output +'.mat'
   output_obj_u = output +'_u.pkl'
   output_obj_v = output +'_v.pkl'
# From here on, always use T,Y,X order
   X = np.concatenate([to,yo,xo],axis=1)
   obs = np.concatenate([vo,uo],axis=1)
   Xt = np.concatenate([tt,yt,xt],axis=1)
   obst = np.concatenate([vt,ut],axis=1)
#########
# Compute covariances
# Pay attention on the order T,Y,X
   priors = np.array([2.,2.,2.])
   bounds = np.array([[0.5,10.],[0.3,10],[0.3,10]])
   if gpy==1:
      kv = GPy.kern.RBF(input_dim=3,ARD=True) #+ GPy.kern.RBF(input_dim=3,ARD=True)
      model = GPy.models.GPRegression(X,vo,kv)
      model.optimize_restarts(messages=False,num_restarts=num_restarts)
      print model
      print kv.parameters
      del model,kv
      print gc.collect(2)

      kv = GPy.kern.RBF(input_dim=3,ARD=True) #+ GPy.kern.RBF(input_dim=3,ARD=True)
      model = GPy.models.GPRegression(X,uo,kv)
      model.optimize_restarts(messages=False,num_restarts=num_restarts)
      print model
      print kv.parameters

   else:
      noise = kernels.WhiteKernel(noise_level=0.001, noise_level_bounds=(1e-06, 1.0))
      kv = kernels.RBF(length_scale=priors) #,length_scale_bounds=bounds)
      k = kernels.Sum(kv,noise)
      model_v = GaussianProcessRegressor(kernel=k,n_restarts_optimizer=num_restarts)
      model_v.fit(X,vo)
      print model_v.kernel_
      del model_v,kv,k
      print gc.collect(2)

      ku = kernels.RBF(length_scale=priors) # ,length_scale_bounds=bounds)
      k = kernels.Sum(ku,noise)
      model_u = GaussianProcessRegressor(kernel=k,n_restarts_optimizer=num_restarts)
      model_u.fit(X,uo)
      print model_u.kernel_

#   model.pickle(output_obj_v)
#   del model,kv
#   print gc.collect(2)

## Compute U
#   ku =  GPy.kern.RBF(input_dim=3,ARD=True) + GPy.kern.RBF(input_dim=3,ARD=True) #k.copy()
#   model = GPy.models.GPRegression(X,uo,ku)
#   model.optimize_restarts(messages=False,num_restarts=num_restarts)
#   model.pickle(output_obj_u)
#   print model
#   print ku.parameters
#   sio.savemat(output_mat,{'Xo':X,'obs':obs,'Xt':Xt,'test_points':obst})
   print 'End of script, time : ' + str(datetime.now()-startTime)
#   print gc.collect(2)

mem_usage = mprof.memory_usage(kriging, timeout=60, interval=0.2)
print mem_usage
print np.max(mem_usage)


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
sys.path.append('/nethome/swang/python/lib/python2.7/site-packages')
from netCDF4 import Dataset

#import GP_scripts as gps

def getData(st,et):
   with open('Filtered_2016_2_7.pkl','rb') as input:
       tr = pickle.load(input)
   # drifter L_0937 has a strage velocities on its first 2-3 steps
   # this drifter is index 238 on file 'Filtered_2016_2_7.pkl'
   # to check, figure(); quiver(tr.lon[8:20,238],tr.lat[8:20,238],tr.u[8:20,238],tr.v[8:20,238])
   # exclude data from this drifter using NaN
   tr.lat[:,238] = np.nan
   tr.lon[:,238] = np.nan
   tr.u[:,238] = np.nan
   tr.v[:,238] = np.nan

   time = (tr.time[st:et] - tr.time[st])/3600. # time in hours    
   latt = tr.lat[st:et,:]
   lont = tr.lon[st:et,:] 
   uob = tr.u[st:et,:]
   vob = tr.v[st:et,:]
   N = np.size(latt,1)
   validPoints = np.zeros(N)
   for i in range(N):
      valid = np.where((~np.isnan(lont[:,i]))&(~np.isnan(latt[:,i])))[0]
      validPoints[i] = np.size(valid)
   
   order = np.squeeze(validPoints.argsort(axis=0))  
   order = order[::-1] # reverse order
   
   return time,latt[:,order],lont[:,order],vob[:,order],uob[:,order],validPoints[order]

########################################################################################
def kriging(st,et,sample_step=5,skip=5,num_restarts=10,save=0,output='rbfModel'):
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
   rx = 1.
   sigx = 0.1
   ry = 1.
   sigy = 0.1
   rt = 2.
   sigt = 0.5
   noise = 0.0002
#########
# Compute covariances
# Pay attention on the order T,Y,X

#
   kt = GPy.kern.RBF(input_dim=1, active_dims=[0], variance=sigt, lengthscale=rt)
   ky = GPy.kern.RBF(input_dim=1, active_dims=[1], variance=sigy, lengthscale=ry)
   kx = GPy.kern.RBF(input_dim=1, active_dims=[2], variance=sigx, lengthscale=rx)
   k = kt * ky * kx 
   print k
   kv = k.copy() 
## Compute V
   model_v = GPy.models.GPRegression(X,vo,kv)
#   model_v.Gaussian_noise = noise # got from previous experiments
#   model_v.optimize()   
#   print model_v
   model_v.optimize_restarts(messages=True,num_restarts=num_restarts)
   print model_v
#   vg,vgVar = model_v.predict(Xrg)
   vHP = model_v.param_array

## Compute U
   ku = k.copy()
   model_u = GPy.models.GPRegression(X,uo,ku)
#   model_u.Gaussian_noise = noise # got from previous experiments
#   model_u.optimize()   
#   print model_u
   model_u.optimize_restarts(messages=True,num_restarts=num_restarts)
   print model_u
#   ug,ugVar = model_u.predict(Xrg)
#   uHP = model_u.param_array
#   uHP_names = model_u.parameter_names()

#   vg = np.reshape(vg,[tg.size,yg.size,-1])
#   ug = np.reshape(ug,[tg.size,yg.size,-1])
#   vgVar = np.reshape(vgVar,[tg.size,yg.size,-1])
#   ugVar = np.reshape(ugVar,[tg.size,yg.size,-1])
   if save==0:
      return X,obs,Xt,obst,model_v,model_u
   else:
      sio.savemat(output_mat,{'Xo':X,'obs':obs,'Xt':Xt,'test_points':obst})
#      with open(output_obj_u,'wb') as output:
#           pickle.dump(model_u,output,-1)
#      with open(output_obj_v,'wb') as output:
#           pickle.dump(model_v,output,-1)
   model_v.pickle(model_obj_v)
   model_u.pickle(model_obj_u)

   print 'End of script, time : ' + str(datetime.now()-startTime)

#######################################################################################
def predict(filename,tlim=[0,0],ylim=[0,0],xlim=[0,0],dt=0.5,dx=0.5,xL=40,yL=40):
   startTime = datetime.now()
   cheatPickle = GPy.load('cheatPickle.pkl')# - Stirr pickle to the right class,
                                            #       otherwise it tries to load a paramz object.
                                            # - Pickle is not working properly with GPy objects 
                                            #       saved on Pegasus 2, it only works with 
                                            #       the ones saved in my WS. 
   model_v = GPy.load(filename + '_v.pkl')
   model_u = GPy.load(filename + '_u.pkl')
#   with open(filename+'_v.pkl','rb') as input:
#        model_v = pickle.load(input)
#   with open(filename+'_u.pkl','rb') as input:
#        model_u = pickle.load(input)
   if (ylim[0] == ylim[1])&(xlim[0] == xlim[1]):
      f = sio.loadmat(filename+'.mat')
      Xo = f['Xo']
      Xp,tp,yp,xp = getGrid(Xo[:,0],Xo[:,1],Xo[:,2])
   else:
      Xp,tp,yp,xp = getGrid(tlim,ylim,xlim,dt,dx,xL,yL)
   inc = yp.size*xp.size
   i2=0
   for i in range(tp.size): 
      Xp2 = Xp[i2:i2+inc,:]
      V2,VVar2 = model_v.predict(Xp2)
      U2,UVar2 = model_u.predict(Xp2)
      if i==0:
         V = V2
         VVar = VVar2
         U = U2
         UVar = UVar2
      else:
         V = np.concatenate([V,V2],axis=0)
         U = np.concatenate([U,U2],axis=0)
         VVar = np.concatenate([VVar,VVar2],axis=0)
         UVar = np.concatenate([UVar,UVar2],axis=0)
      i2+=inc
      print 'step '+str(i+1)+'/'+str(tp.size)
#   V,VVar = model_v.predict(Xp)
#   U,UVar = model_u.predict(Xp)
   print 'Creating Netcdf file; running time = ' + str(datetime.now()-startTime)
   createNC(filename+'_2.nc',tp,yp,xp)
   V = np.reshape(V,[tp.size,yp.size,xp.size])
   VVar = np.reshape(VVar,[tp.size,yp.size,xp.size])
   U = np.reshape(U,[tp.size,yp.size,xp.size])
   UVar = np.reshape(UVar,[tp.size,yp.size,xp.size])
   fi = Dataset(filename+'_2.nc','a')
   fi = writeNC(fi,'v',V)
   fi = writeNC(fi,'u',U)
   fi = writeNC(fi,'vvar',VVar)
   fi = writeNC(fi,'uvar',UVar)
   fi.close()
   print 'End of script, time : ' + str(datetime.now()-startTime)

    #return Xp,V,U,VVar,UVar
#######################################################################################
def predictTest(filename):
   startTime = datetime.now()
   cheatPickle = GPy.load('cheatPickle.pkl')# - Stirr pickle to the right class,
                                            #       otherwise it tries to load a paramz object.
                                            # - Pickle is not working properly with GPy objects 
                                            #       saved on Pegasus 2, it only works with 
                                            #       the ones saved in my WS. 
   model_v = GPy.load(filename + '_v.pkl')
   model_u = GPy.load(filename + '_u.pkl')
   f = sio.loadmat(filename+'.mat')
   Xt = f['Xt']
   obst=f['test_points']
   vt = obst[:,0]
   ut = obst[:,1]
   Nt = ut.size()
   step = Nt/10
   
   for i in range(0,Nt,step):
      if i+step<Nt:
         Xt2 = Xt[i:i+step,:]
      elif i<Nt:
         Xt2 = Xt[i:,:]
      V2,VVar2 = model_v.predict(Xt2)
      U2,UVar2 = model_u.predict(Xt2)
      if i==0:
         V = V2
         VVar = VVar2
         U = U2
         UVar = UVar2
      else:
         V = np.concatenate([V,V2],axis=0)
         U = np.concatenate([U,U2],axis=0)
         VVar = np.concatenate([VVar,VVar2],axis=0)
         UVar = np.concatenate([UVar,UVar2],axis=0)
      print 'step '+str(i+1)+'/'+str(10)

   print 'End of script, time : ' + str(datetime.now()-startTime)
   output = filename +'_test.mat'
   sio.savemat(output,{'Xt':Xt,'Vp':V,'VpVar':VVar,'Up':U,'UpVar':UVar,'test_points':obst})

    #return Xp,V,U,VVar,UVar

#######################################################################################
    
def getRMSE(filename):
    with open(filename+'_v.pkl','rb') as input:
         model_v = pickle.load(input)
    with open(filename+'_u.pkl','rb') as input:
         model_u = pickle.load(input)

    f = sio.loadmat(filename+'.mat')
    Xt = f['Xt']
    obst = f['test_points']
    vt = obst[:,0]
    ut = obst[:,1]
    rmse_v = testModel1D(model_v,Xt,vt)
    print rmse_v
    rmse_u = testModel1D(model_u,Xt,ut)
    print rmse_u
    return rmse_v,rmse_u
###########################################################################################
def testModel1D(model,Xt,test_points):
    f,fVar = model.predict(Xt)
    return rmse(f,test_points)
#################################################################################
def rmse(ys,y):
   # compute the rmse 
   error = ys-y
   error = np.reshape(error,[error.size])
   return np.sqrt(np.mean(np.square(error)))

######################################################################################
def getGrid(to,yo,xo,dt=0.5,dx=0.5,xL=40,yL=40):
# GRID check size of the final matrix Xrg (reshaped grid)
   if (xo.max()-xo.min())>xL:
      xmin = xo.mean() - xL/2
      xmax = xo.mean() + xL/2
      print 'here x'
   else:
      xmin = xo.min()- dx
      xmax = xo.max() + dx
 
   if (yo.max()-yo.min())>yL:
      ymin = yo.mean() - yL/2
      ymax = yo.mean() + yL/2
      print 'here y'
   else:
      ymin = yo.min() - dx
      ymax = yo.max() + dx

   xg = np.arange(xmin,xmax,dx)
   yg = np.arange(ymin,ymax,dx)
   tg = np.arange(to.min(),to.max(),dt)

   Yg,Tg,Xg = np.meshgrid(yg,tg,xg) # Works
   #  1st index vary with T
   #  2nd index vary with Y
   #  3rd index vary with X

   Tr = np.reshape(Tg,[Tg.size,1])
   Yr = np.reshape(Yg,[Yg.size,1]) 
   Xr = np.reshape(Xg,[Xg.size,1])
   return np.concatenate([Tr,Yr,Xr],axis=1),tg,yg,xg





import numpy as np
from datetime import datetime,timedelta
import cPickle as pickle
import GPy
from sklearn.gaussian_process import kernels,GaussianProcessRegressor
import os
import pyproj
from matplotlib import rc,rcParams
import myKernel2
import scipy.io as sio
from printNCFiles import createNC,writeNC
import sys
from netCDF4 import Dataset
import gc
#import GP_scripts as gps
#from sklearn.gaussian_process import kernels,GaussianProcessRegressor
lat0 = 28.8
lon0 = -88.6
NAD83=pyproj.Proj("+init=EPSG:3452") #Louisiana South (ftUS)
x_ori,y_ori=NAD83(lon0,lat0)
################################################################
def saveDict(obj, filename):
    if np.size(obj) == 1:
       with open(filename, 'wb') as output:
           pickle.dump(obj, output)
    else:
       with open(filename, 'wb') as output:
           for i in range(np.size(obj)):
               pickle.dump(obj[i], output)
#################################################################
def read_object(filename):
    lists = []
    infile = open(filename, 'r')
    while True:
        try:
            lists.append(pickle.load(infile))
        except (EOFError, pickle.UnpicklingError):
            break
    infile.close()
    return lists
#####################################################################################################
def getData(st,et,laser=1):
   if laser==1:
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
      time = (tr.time[st:et] - tr.time[0])/3600. # time in hours    
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
   else:
      f = Dataset('Simulations/Output.nc','r')
      print st,et
      latt = f.variables['lat'][st:et][:]
      lont = f.variables['lon'][st:et][:]
      vob = f.variables['v'][st:et][:]
      uob = f.variables['u'][st:et][:]
      initday = f.variables['time'][0]
      time=  f.variables['time'][st:et]-initday
      validPoints = np.zeros(np.size(latt,1)) + time.size
      order = range(np.size(latt,1)) # don't do anything here
   return time,latt[:,order],lont[:,order],vob[:,order],uob[:,order],validPoints[order]
############################################################################################
def boundData(var,varlim,lat,lon,v,u):
    """
    Bound initial position based on var value.
    Var can be lat, lon, x, y, u and v
    Variables shoul be organaized as 2D arrays (time X number of drifters)
    """
    il = np.where((var[0,:]>=varlim[0])&(var[0,:]<=varlim[1]))[0]
    return lat[:,il], lon[:,il],v[:,il], u[:,il]
############################################################################################
def scikit_prior(filename0,varname='v',dt=0,tlim=6,radar='',xlim=[0,0],ylim=[0,0],dx=0,ind=0,xrange=3):
   startTime = datetime.now()
   dir0,a = filename0.split("res")
   b,fname0 = a.split("/") 
   fname0 = dir0 + fname0
   fm = sio.loadmat(fname0 + '.mat') 
   print 'Longitude limits:',  xlim
   print 'Latitude limits :',  ylim

   # get radar data grid, if that is the case:
   if radar!='':
      inFile = Dataset(radar, 'r')
      lon0,lat0 = inFile.variables['imageOriginPosition'][:]
      x0,y0=NAD83(lon0,lat0)
      x0 = (x0 - x_ori)/1000. # in km
      y0 = (y0 - y_ori)/1000. # in km
      xg = x0 + inFile.variables['xCoords'][:]/1000.
      yg = y0 + inFile.variables['yCoords'][:]/1000.
      tr = inFile.variables['time'][:]
      ur = inFile.variables['ux'][:]
      vr = inFile.variables['uy'][:]
      t0 = datetime(2016,01,01) # radar data counts from here, in hours
      t0D = datetime(2016, 2, 7, 2, 15) # first time from Filtered_2016_2_7.pkl'
      tg = np.array([(t0 + timedelta(tr[0]) - t0D).total_seconds()/3600])
      it=0
      filename = filename0 + '_radar'
      Yg,Tg,Xg = np.meshgrid(yg,tg,xg)
      Tg = np.reshape(Tg,[Tg.size,1])
      Yg = np.reshape(Yg,[Yg.size,1]) 
      Xg = np.reshape(Xg,[Xg.size,1])
      X = np.concatenate([Tg,Yg,Xg],axis=1)
   else:   # DEFINE GRID
      if (xlim[1]>xlim[0])&(ylim[1]>ylim[0]): # should focus here
         X,tcenter,yg,xg = getGrid([dt,dt+1],ylim,xlim,1,dx)
         filename = filename0 + '_cyc'
      else: # this is for preexisting grids
         f = Dataset(filename0 + '.nc','r')
#         HPU = f.variables['hyperparam_u'][:]
#         HPV = f.variables['hyperparam_v'][:]
         xg = f.variables['x'][:]
         yg = f.variables['y'][:]
         tg = f.variables['time'][:]
         it = dt #tg.size/2 + dt
         tcenter = np.array([tg[it]]) 
         Yg,Tg,Xg = np.meshgrid(yg,tg,xg)
         Tg = np.reshape(Tg,[Tg.size,1])
         Yg = np.reshape(Yg,[Yg.size,1]) 
         Xg = np.reshape(Xg,[Xg.size,1])
         X = np.concatenate([Tg,Yg,Xg],axis=1)
         filename = filename0
         inc = yg.size * xg.size 
         i2= inc*it
         X = X[i2:i2+inc,:]

   filename = filename +'_'+ str(np.round(tcenter[0],decimals=2))  + 'h_scikit_'
   outFile = filename +str(ind)+'.nc'

# LOAD Observations  
   to = fm['Xo'][:,0]   
   tt = fm['Xt'][:,0]   
   xo = fm['Xo'][:,2]   
   xt = fm['Xt'][:,2]   
   ito = np.where((to>=tcenter-tlim)&(to<=tcenter+tlim)&(xo>=xlim[0]-xrange)&(xo<=xlim[1]+xrange))
   itt = np.where((tt>=tcenter-tlim)&(tt<=tcenter+tlim)&(xt>=xlim[0]-xrange)&(xt<=xlim[1]+xrange))
   Xo = fm['Xo'][ito,:].squeeze()
   Xt = fm['Xt'][itt,:].squeeze()
   XT = np.concatenate([Xo,Xt],axis=0)
   print 'Number of observation points: ',np.size(XT,0)
   obs = fm['obs'][ito,:].squeeze()   
   obst = fm['test_points'][itt,:].squeeze()   
# LOAD Hyper-Parameters
   cheatPickle = GPy.load('cheatPickle.pkl')
   model = GPy.load(fname0 +'_'+varname+'.pkl')
   HP = model.param_array
   covarname = varname + 'var'
   modelName = filename + varname + '.pkl'
   if varname=='u':
      u = np.concatenate([obs[:,1],obst[:,1]])[:,None]
   else:
      u = np.concatenate([obs[:,0],obst[:,0]])[:,None]
   N = HP.size - 1
   noise = HP[-1]
   print 'noise = ' + str(HP[-1])
# Build Model
   print modelName   
#   if not os.path.isfile(modelName):
   k = HP[0]* kernels.RBF(length_scale=[HP[1],HP[2],HP[3]])
   print 'var1 = '+str(HP[0])
   if N > 5:
      i=4
      k = k + HP[i]* kernels.RBF(length_scale=[HP[i+1],HP[i+2],HP[i+3]])
      print 'var2 = ' + str(HP[i])    
   k  = k + kernels.WhiteKernel(noise_level=noise)
   print k
   model_u = GaussianProcessRegressor(kernel=k,optimizer=None)
   print np.size(XT,0),np.size(XT,1)
   print np.size(u,0), np.size(u,1)
   model_u.fit(XT,u)
# file might be to large to save
#      with open(modelName,'wb') as output:
#      pickle.dump(model_u,open(modelName,'wb'))
#   else:
#      with open(modelName,'rb') as input:
#           model_u = pickle.load(input)
 
# REGRESSION   
   U,Ustd = model_u.predict(X,return_std=True)
   U = np.reshape(U,[tcenter.size,yg.size,xg.size])
   Ustd = np.reshape(Ustd,[tcenter.size,yg.size,xg.size])
# SAVE NETCDF
   if not os.path.isfile(outFile):
      createNC(outFile,tcenter,yg,xg,HP)
   print np.ndim(U),np.size(U,0),np.size(U,1)
   print np.ndim(Ustd),np.size(Ustd,0),np.size(Ustd,1)
   fi = Dataset(outFile,'a')
   fi = writeNC(fi,varname,U) 
   fi = writeNC(fi,covarname,Ustd**2)
   fi = writeNC(fi,'hyperparam_'+varname,HP)
   fi.close() 
   print 'End of script, time : ' + str(datetime.now()-startTime)
        
############################################################################################
def scikitSnapshot(filename,var,dt,ylim=[0,0],xlim=[0,0]):
   f = Dataset(filename + '.nc','r')
   HPU = f.variables['hyperparam_u'][:]
   HPV = f.variables['hyperparam_v'][:]
   xg = f.variables['x'][:]
   yg = f.variables['y'][:]
   tg = f.variables['time'][:]

   if (xlim[1]>xlim[0])&(ylim[1>ylim[0]]):
      xg = xg[np.where((xg>=xlim[0])&(xg<=xlim[1]))]
      yg = yg[np.where((yg>=ylim[0])&(yg<=ylim[1]))]

   Yg,Tg,Xg = np.meshgrid(yg,tg,xg)
   Tg = np.reshape(Tg,[Tg.size,1])
   Yg = np.reshape(Yg,[Yg.size,1]) 
   Xg = np.reshape(Xg,[Xg.size,1])
   X = np.concatenate([Tg,Yg,Xg],axis=1)
   inc = yg.size * xg.size 

   it = tg.size/2 + dt
   i2= inc*it
   outFile = filename +'_'+ str(tg[it])  + 'h_scikit.nc'
   tg = np.array([tg[it]])

   X2 = X[i2:i2+inc,:]
   print np.size(X2,0),np.size(X2,1)

   if var == 1: 
      with open(filename + '_scikit_u.pkl','rb') as input:
           model_u = pickle.load(input) 
   else:
      with open(filename + '_scikit_v.pkl','rb') as input:
           model_u = pickle.load(input)

   U,Ustd = model_u.predict(X2,return_std=True)

   print tg.size
   U = np.reshape(U,[tg.size,yg.size,xg.size])
   Ustd = np.reshape(Ustd,[tg.size,yg.size,xg.size])

   if not os.path.isfile(outFile):
      createNC(outFile,tg,yg,xg,HPV)
   fi = Dataset(outFile,'a')
   fi = writeNC(fi,varname,U) 
   fi = writeNC(fi,covarname,Ustd**2)
   fi.close() 
   print 'End of script, time : ' + str(datetime.now()-startTime)

#########################################################################################
def kriging(st,et,lalim=[0,0],lolim=[0,0],sample_step=5,skip=5,nKernels=1,output='rbfModel',pkg='GPy',kernelType = 1,laser=1):
   """
   Estimate velocity field using rbf kernels on t,y,x (K = Kt*Ky*Kx)
   st,et = initial,final time step of drifter data
   lalim,lolim = lat lon limits for initial data. If none, use all data
   sample_step = pick data each sample_step. If negative, pick data at 
each time step = sample_step*(-1)     
   skip = skip drifters. If skip ==1, use all drifters
   nKernels = the kernel will be a sum of nKernels of the same type.
   kernelType = 1 for rbf, 2 for divergence-free rbf, 3 for curl-free 
rbf or 4 for divergence-free+curl-free rbf.
   laser = 1 to use laser data. laser!=1 will use numerical model outputs
   """
#   st = 0; et = 144
   startTime = datetime.now()
   time,latt,lont,vob,uob,validPoints = getData(st,et,laser)

   # Next is to select drifters that at timestep st are bounded by lolim and/or lalim
   # This part was designed to select data around the radar velocity estimates 
   if (lolim[1]>lolim[0]):
      latt,lont,vob,uob =  boundData(lont,lolim,latt,lont,vob,uob)
      print 'data limited by initial longitude, between ' + str(lolim[0])+' and '+str(lolim[1])
      print 'Total number of data points selected: '+ str((np.size(latt))) 
   if (lalim[1]>lalim[0]):
      latt,lont,vob,uob =  boundData(latt,lalim,latt,lont,vob,uob)
      print 'data limited by initial latitude, between ' + str(lalim[0])+' and '+str(lalim[1])
      print 'Total number of data points selected: '+ str((np.size(latt))) 

   # origin of cartesian coord.
#   lat0 = 28.8
#   lon0 = -88.6
#   NAD83=pyproj.Proj("+init=EPSG:3452") #Louisiana South (ftUS)
   xob,yob=NAD83(lont,latt)
   tob = time[:,None]; tob = np.repeat(tob,np.size(latt,1),axis=1)
   yob[np.where(np.isnan(lont))]=np.nan
   xob[np.where(np.isnan(lont))]=np.nan

#   x_ori = np.nanmin(xob)+2; y_ori = np.nanmin(yob)+2
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
      lat_o = np.reshape(latt[samples,::skip],[-1,1])
      lon_o = np.reshape(lont[samples,::skip],[-1,1])
      uo = np.reshape(uob[samples,::skip],[-1,1])
      vo = np.reshape(vob[samples,::skip],[-1,1])

      tt    = np.reshape(tob[testt[:,None],testd],[-1,1])
      yt    = np.reshape(yob[testt[:,None],testd],[-1,1])
      xt    = np.reshape(xob[testt[:,None],testd],[-1,1])
      lat_t = np.reshape(latt[testt[:,None],testd],[-1,1])
      lon_t = np.reshape(lont[testt[:,None],testd],[-1,1])
      ut    = np.reshape(uob[testt[:,None],testd],[-1,1])
      vt    = np.reshape(vob[testt[:,None],testd],[-1,1])

   else:
      ss = 0
      samples = np.arange(0,xob.size,sample_step) #np.random.randint(0,xo.size,nsamples)
      test = set(np.arange(xob.size)) - set(samples)
      test = np.array(list(test))
      xt = np.reshape(xob,[-1])[test,None]
      yt = np.reshape(yob,[-1])[test,None]
      lat_t = np.reshape(latt,[-1])[test,None]
      lon_t = np.reshape(lont,[-1])[test,None]
      tt = np.reshape(tob,[-1])[test,None]
      ut = np.reshape(uob,[-1])[test,None]
      vt = np.reshape(vob,[-1])[test,None]
      xo = np.reshape(xob,[-1])[samples,None]
      yo = np.reshape(yob,[-1])[samples,None]
      lon_o = np.reshape(lont,[-1])[samples,None]
      lat_o = np.reshape(latt,[-1])[samples,None]
      to = np.reshape(tob,[-1])[samples,None]
      uo = np.reshape(uob,[-1])[samples,None] 
      vo = np.reshape(vob,[-1])[samples,None]

  
   validPoints = np.where((~np.isnan(xo))&(~np.isnan(yo)))
   to = to[validPoints][:,None]
   xo = xo[validPoints][:,None]
   yo = yo[validPoints][:,None]
   lon_o = lon_o[validPoints][:,None]
   lat_o = lat_o[validPoints][:,None]
   uo = uo[validPoints][:,None]
   vo = vo[validPoints][:,None]
   validPoints = np.where((~np.isnan(xt))&(~np.isnan(yt)))
   tt = tt[validPoints][:,None]
   xt = xt[validPoints][:,None]
   yt = yt[validPoints][:,None]
   lon_t = lon_t[validPoints][:,None]
   lat_t = lat_t[validPoints][:,None]
   ut = ut[validPoints][:,None]
   vt = vt[validPoints][:,None]

   print 'number of observations: '+str(np.size(vo))
   output_mat = output +'.mat'
   output_obj_u = output +'_u.pkl'
   output_obj_v = output +'_v.pkl'
# From here on, always use T,Y,X order
   X = np.concatenate([to,yo,xo],axis=1)
   LL_o = np.concatenate([to,lat_o,lon_o],axis=1)
   obs = np.concatenate([vo,uo],axis=1)
   Xt = np.concatenate([tt,yt,xt],axis=1)
   LL_t = np.concatenate([tt,lat_t,lon_t],axis=1)
   obst = np.concatenate([vt,ut],axis=1)
#########
# Compute covariances
# Pay attention on the order T,Y,X
#
#   if pkg=='GPy': #use Gpy package
   if kernelType==1:
      k2 = GPy.kern.RBF(input_dim=3,ARD=True) 
      Obs = [vo,uo,':P']
      output_obj = [output_obj_v,output_obj_u]
   else:
      obs = np.concatenate([vo,uo],axis=0)
      obst = np.concatenate([vt,ut],axis=0)
      Obs = [obs,':P']
    
      if kernelType==2:
         k2 = myKernel2.divFreeK(input_dim=3, active_dims=[0,1,2], var=1., lt=1., ly=1., lx=1.)
         output_obj = [output+'_divFree.pkl']
      elif kernelType==3:
         output_obj = [output+'_curlFree.pkl']
         k2 = myKernel2.curlFreeK(input_dim=3, active_dims=[0,1,2], var=1., lt=1., ly=1., lx=1.)
      else:
         output_obj = [output+'_combined.pkl']
         k2 = myKernel2.divFreeK(input_dim=3) + myKernel2.curlFreeK(input_dim=3)
   k = k2.copy()
   for i in range(nKernels-1):
       k = k + k2
   print k
   for i in range(np.size(Obs)-1):
       kv = k.copy()
       model = GPy.models.GPRegression(X,Obs[i],k)
       model.pickle(output_obj[i]) 
       del model,kv
       print gc.collect(2) # delete garbage


   sio.savemat(output_mat,{'Xo':X,'obs':obs,'Xt':Xt,'LL_o':LL_o,'LL_t':LL_t,'test_points':obst})
   print 'End of script, time : ' + str(datetime.now()-startTime)

#   else:
#      model_v = GaussianProcess(theta0=[1,1,1],thetaL=[0.1,0.1,0.1],thetaU=[10,10,10],nugget = 0.0001,random_start=10)
#      model_v.fit(X,vo)
#      with open('scikit_' + output_obj_v,'wb') as output:
#           pickle.dump(model_v,output,-1)
#      model_u = GaussianProcess(theta0=[1,1,1],thetaL=[0.1,0.1,0.1],thetaU=[10,10,10],nugget = 0.0001,random_start=10)
#      model_u.fit(X,uo)
#      with open('scikit_' + output_obj_u,'wb') as output:
#           pickle.dump(model_u,output,-1)
#######################################################################################
def runRestarts(fname,nres=10,nKernels=2):
   filename = fname+'.pkl'
   startTime = datetime.now()
   cheatPickle = GPy.load('cheatPickle.pkl')# - Stirr pickle to the right class,
                                            #       otherwise it tries to load a paramz object.
                                            # - Pickle is not working properly with GPy objects 
                                            #       saved on Pegasus 2, it only works with 
                                            #       the ones saved in my WS. 
   model = GPy.load(filename) 
   if (len(model.optimization_runs)>0):
   # This part is to overcome a bug when reopening a model that was  previously 
   # optimized. For some reason, the dictionaries of the optimization_runs are lost 
   # when the model is pickled. 
   # So I'm saving these dictionaries in a list as 'dict_'+filename, and loading 
   # them whenever I want to carry out more optimization runs. 
      optDict = read_object(fname + '_dict.pkl')
      for i in range(len(model.optimization_runs)):
          model.optimization_runs[i].__dict__ = optDict[i]
   ###
   hyp_old = model.param_array[:]
   model.optimize_restarts(messages=False,num_restarts=nres)
   hyp = model.param_array
   model.pickle(filename)
   optDict = []
   ### saving dictionaries of the optimization runs
   for i in range(len(model.optimization_runs)):
       optDict.append(model.optimization_runs[i].__dict__)
   saveDict(optDict,fname+'_dict.pkl')
   print 'Optimized Hyperparameters ================================================='
   for i in range(nKernels):
       print 'Var.' + str(i+1) + ' = ' + str(hyp_old[4*i]) + '   :   ' + str(hyp[4*i]) 
       print 'Lt ' + str(i+1) + '  = ' + str(hyp_old[4*i+1]) + '   :   ' + str(hyp[4*i+1])
       print 'Ly ' + str(i+1) + '  = ' + str(hyp_old[4*i+2]) + '   :   ' + str(hyp[4*i+2]) 
       print 'Lx ' + str(i+1) + '  = ' + str(hyp_old[4*i+3]) + '   :   ' + str(hyp[4*i+3])

   print 'Noise = ' + str(hyp_old[-1]) + '   :   ' + str(hyp[-1]) 
   print '==========================================================================='

   print 'End of script, time : ' + str(datetime.now()-startTime)

#######################################################################################
def predict(filename,tlim=[0,0],ylim=[0,0],xlim=[0,0],dt=0.5,dx=0.5,xL=40,yL=40,Simul=0):
   startTime = datetime.now()
   cheatPickle = GPy.load('cheatPickle.pkl')# - Stirr pickle to the right class,
                                            #       otherwise it tries to load a paramz object.
                                            # - Pickle is not working properly with GPy objects 
                                            #       saved on Pegasus 2, it only works with 
                                            #       the ones saved in my WS. 
   model_v = GPy.load(filename + '_v.pkl')
   hypv = model_v.param_array
   hypv_names = model_v.parameter_names()

   model_u = GPy.load(filename + '_u.pkl')
   hypu = model_u.param_array
#   hypu_names = model_u.parameter_names()

#   with open(filename+'_v.pkl','rb') as input:
#        model_v = pickle.load(input)
#   with open(filename+'_u.pkl','rb') as input:
#        model_u = pickle.load(input)
      

   if (ylim[0] == ylim[1])&(xlim[0] == xlim[1]):
      f = sio.loadmat(filename+'.mat')
      Xo = f['Xo']
      print Xo[:,0].min() ,Xo[:,0].max() 
      print Xo[:,1].min() ,Xo[:,1].max() 
      print Xo[:,2].min() ,Xo[:,2].max() 

      if Simul==1: # get NCOM grid
         ncom = Dataset( 'osprein_2013_8.nc', 'r')
         ylim = np.array([Xo[:,1].min() ,Xo[:,1].max()])
         xlim = np.array([Xo[:,2].min() ,Xo[:,2].max()])
         lon_nc = ncom.variables['lon'][:]
         lat_nc = ncom.variables['lat'][:]
         print lon_nc.size,lat_nc.size
         ti_nc = ncom.variables['time'][1:]-ncom.variables['time'][1]
         dt = ti_nc[1]
         tp = np.arange(Xo[:,0].min(),Xo[:,0].max()+dt,dt)
         Lon_nc,Lat_nc = np.meshgrid(lon_nc,lat_nc)
         xnc,ync=NAD83(Lon_nc,Lat_nc)
         xnc2 = (np.reshape(xnc,[1,-1])-x_ori)/1000.
         ync2 = (np.reshape(ync,[1,-1])-y_ori)/1000.
         limit = np.where((xnc2>=xlim[0])&(xnc2<=xlim[1])&(ync2>=ylim[0])&(ync2<=ylim[1]))
         # The next part is necessary for the netcdf grid, so that the space dimensions 
         # can be defined by the model grid's lon and lat 
         Lo_mx = np.reshape(Lon_nc,[1,-1])[limit].max() # get the lon inside the space limit that we want
         Lo_mn = np.reshape(Lon_nc,[1,-1])[limit].min() # get the lon inside the space limit that we want
         La_mx = np.reshape(Lat_nc,[1,-1])[limit].max()
         La_mn = np.reshape(Lat_nc,[1,-1])[limit].min()
         # make the min,max of lon and lat to define the limits of the grid, and not x and y
         limit2 = np.where((Lon_nc>=Lo_mn)&(Lon_nc<=Lo_mx)&(Lat_nc>=La_mn)&(Lat_nc<=La_mx))
         xnc = (xnc[limit2][None,:]-x_ori)/1000. 
         xnc = np.repeat(xnc,tp.size,axis=0)
         xnc = np.reshape(xnc,[-1,1])
         ync = (ync[limit2][None,:]-y_ori)/1000. 
         ync = np.repeat(ync,tp.size,axis=0)
         ync = np.reshape(ync,[-1,1])
         tnc = np.repeat(tp[:,None],xnc.size/tp.size,axis=1)         
         tnc = np.reshape(tnc,[-1,1])
         Xp = np.concatenate([tnc,ync,xnc],axis=1)
         # xp,yp are lon,lat here!
         xp = lon_nc[np.where((lon_nc>=Lo_mn)&(lon_nc<=Lo_mx))]
         yp = lat_nc[np.where((lat_nc>=La_mn)&(lat_nc<=La_mx))]
         sio.savemat(filename+'_grid.mat',{'Xp':Xp})
      else:
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
      print gc.collect(2)
      print 'step '+str(i+1)+'/'+str(tp.size)
#   V,VVar = model_v.predict(Xp)
#   U,UVar = model_u.predict(Xp)
   print 'Creating Netcdf file; running time = ' + str(datetime.now()-startTime)
   createNC(filename+'.nc',tp,yp,xp,hypv)
   V = np.reshape(V,[tp.size,yp.size,xp.size])
   VVar = np.reshape(VVar,[tp.size,yp.size,xp.size])
   U = np.reshape(U,[tp.size,yp.size,xp.size])
   UVar = np.reshape(UVar,[tp.size,yp.size,xp.size])
   fi = Dataset(filename+'.nc','a')
   fi = writeNC(fi,'v',V)
   fi = writeNC(fi,'u',U)
   fi = writeNC(fi,'vvar',VVar)
   fi = writeNC(fi,'uvar',UVar)
   fi = writeNC(fi,'hyperparam_v',hypv)
   fi = writeNC(fi,'hyperparam_u',hypu)
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
   Nt = ut.size
   step = Nt/10
   
   for i in range(0,Nt,step):
      if (i+step<Nt):
         Xt2 = Xt[i:i+step,:]
      elif (i<Nt):
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
   if (np.max(xo)-np.min(xo))>xL:
      xmin = np.mean(xo) - xL/2
      xmax = np.mean(xo) + xL/2
      print 'here x'
   else:
      xmin = np.min(xo) - dx
      xmax = np.max(xo) + dx
 
   if (np.max(yo)-np.min(yo))>yL:
      ymin = np.mean(yo) - yL/2
      ymax = np.mean(yo) + yL/2
      print 'here y'
   else:
      ymin = np.min(yo) - dx
      ymax = np.max(yo) + dx

   xg = np.arange(xmin,xmax,dx)
   yg = np.arange(ymin,ymax,dx)
   tg = np.arange(np.min(to),np.max(to),dt)

   Yg,Tg,Xg = np.meshgrid(yg,tg,xg) # Works
   #  1st index vary with T
   #  2nd index vary with Y
   #  3rd index vary with X

   Tr = np.reshape(Tg,[Tg.size,1])
   Yr = np.reshape(Yg,[Yg.size,1]) 
   Xr = np.reshape(Xg,[Xg.size,1])
   return np.concatenate([Tr,Yr,Xr],axis=1),tg,yg,xg





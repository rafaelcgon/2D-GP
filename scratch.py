import numpy as np
from datetime import datetime
import cPickle as pickle
import GPy
import mpl_toolkits.basemap.pyproj as pyproj
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import rc,rcParams
import myKernel
rc('text',usetex=True)
import scipy.io as sio
from laser_io_methods import save_object
import GP_scripts as gps
# document to modify or create new scripts
# All scripts should be copied to their original .py file
class laser(object):
    def __init__(self,drifters):
        self.id = drifters[0]
        self.datetime = drifters[1] 
        self.time = drifters[2]
        self.lat = drifters[3]
        self.lon = drifters[4]
        self.posErr = drifters[5]
        self.u = drifters[6]
        self.v = drifters[7]
        self.velErr = drifters[8]

def storeDrifter(drif,did,ddatetime,dtime,dlat,dlon,dposErr,du,dv,dvelErr):
   print did
   did = np.array([did])[:,None]
   ddatetime = np.array(ddatetime)#[:,None]
   dtime = np.array(dtime)[:,None]
   dlat = np.array(dlat)[:,None]
   dlon = np.array(dlon)[:,None]
   dposErr = np.array(dposErr)[:,None]
   du = np.array(du)[:,None]
   dv = np.array(dv)[:,None]
   dvelErr = np.array(dvelErr)[:,None]
   if np.size(drif) > 0:
      drif[0]=np.concatenate([drif[0],did],axis=1)   
      drif[1].append(ddatetime)  
      N = np.size(drif[0])
      print N
      # create timestep
      #print drif[2]
      timemin = np.min([np.min(dtime),drif[2][0]])
      timemax = np.max([np.max(dtime),drif[2][-1]])
      dt = drif[2][1]-drif[2][0]
      new_time = np.arange(timemin,timemax+dt,dt)
      new_lat = np.zeros((new_time.size,N)) + np.nan
      new_lon = np.zeros((new_time.size,N)) + np.nan
      new_dposErr = np.zeros((new_time.size,N)) + np.nan
      new_du = np.zeros((new_time.size,N)) + np.nan
      new_dv = np.zeros((new_time.size,N)) + np.nan
      new_dvelErr = np.zeros((new_time.size,N)) + np.nan

      ind = np.where((new_time>=drif[2][0])&(new_time<=drif[2][-1]))
#      print np.size(drif[3],0),np.size(drif[3],1)
#      print np.size(dlat,0),np.size(dlat,1)
      new_lat[ind,0:N-1] = drif[3]
      new_lon[ind,0:N-1] = drif[4]
      new_dposErr[ind,0:N-1] = drif[5]
      new_du[ind,0:N-1] = drif[6]
      new_dv[ind,0:N-1] = drif[7]
      new_dvelErr[ind,0:N-1] = drif[8]
         
      ind2 = np.where((new_time>=dtime[0])&(new_time<=dtime[-1]))
#      print np.size(new_lat[ind2,N-1],0),np.size(new_lat[ind2,N-1],1)
      new_lat[ind2,N-1] = dlat.T
      new_lon[ind2,N-1] = dlon.T
      new_dposErr[ind2,N-1] = dposErr.T
      new_du[ind2,N-1] = du.T
      new_dv[ind2,N-1] = dv.T
      new_dvelErr[ind2,N-1] = dvelErr.T
      drif = [drif[0],drif[1],new_time,new_lat,new_lon,new_dposErr,new_du,new_dv,new_dvelErr]
#      print N,did, np.size(new_time)
   else:
      print 'here'
      drif = [did,[ddatetime],dtime,dlat,dlon,dposErr,du,dv,dvelErr]
      print np.size(did),did, np.size(dtime)
   return drif
##########################################################################  
def initList():
   return [],[],[],[],[],[],[],[]
##########################################################################  
def readFilteredTracks(): # laser_io_methods.py
# Read trajectories from drifters, create objects for each drifter and save all in output
# ID  Date Time   Latitude  Longitude  Pos Error  U (+E-W) V (+N-S) vel Error
# L_0004 2016-01-21 18:45:00.216001  29.03776831107 -87.68717800671    10.9 0.016 0.389 0.033
    initdate = datetime(2016,02,07,0,0,0)
#    output = 'Filtered'+'_'+str(date.year)+'_'+str(date.month)+'_'+str(date.day)+'.pkl'
#    f = open('/home/rgoncalves/LagOil/LASER/Filtered_data/carthe_laser_spot_drifters_clean_v01.dat')
    f = open('carthe_laser_spot_drifters_clean_v01.dat')
    Lines = f.readlines()

    dr_id = ''
    dr_datetime,dr_time,dr_lat,dr_lon,dr_posErr,dr_u,dr_v,dr_velErr = initList()
    valid_drifter = 0
    last_id = ''
#    Lines = f.readlines()
    drifters = []
    n=1
    n2=0
    n3=0
    #print n
    #print Lines[0]
    for line in Lines:
        if (line[0]!='%'):
           n2+=1
           values = [value for value in line.split(' ')]
           if (values[0] != last_id)&(valid_drifter==1):
              print n3,'/',n2,dr_id,values[0]
              drifters = storeDrifter(drifters,dr_id,dr_datetime,
                                      dr_time,dr_lat,dr_lon,dr_posErr,dr_u,dr_v,dr_velErr)
              dr_datetime,dr_time,dr_lat,dr_lon,dr_posErr,dr_u,dr_v,dr_velErr = initList()

           ## get datetime from line 
           dateString = [dateS for dateS in values[1].split('-')]
           yr = np.int(dateString[0])
           mo = np.int(dateString[1])
           da = np.int(dateString[2])
           timeString = [timeS for timeS in values[2].split(':')]
           hr = np.int(timeString[0])
           mi = np.int(timeString[1])
           se = np.int(np.floor(np.float(timeString[2]))) 
           ms = np.int((np.float(timeString[2])-se)*1000000)
           drdate = datetime(yr,mo,da,hr,mi) #,se,ms)
           ## check if new line is from a different drifter
           if ((values[0] != last_id)|(valid_drifter==1))&(values[0][0]=='L')&(drdate >= initdate):
#              print n, values[0], values[3], values[4]
              n3 +=1 
              valid_drifter = 1   
              dr_id = values[0]
              dr_datetime.append(drdate)
              dr_time.append((drdate-initdate).total_seconds())
              dr_lat.append(np.float(values[3]))
              dr_lon.append(np.float(values[4]))
              dr_posErr.append(np.float(values[5]))
              dr_u.append(np.float(values[6]))
              dr_v.append(np.float(values[7]))
              dr_velErr.append(np.float(values[8]))
              if n3==1:
                 print values[0],drdate  
           else: 
              valid_drifter = -1
#              if (values[0][0]!= 'U'):
#                 print values[0],drdate 
           last_id = values[0]
           if (line==Lines[-1])&(valid_drifter==1):
              drifters = storeDrifter(drifters,dr_id,dr_datetime,
                                      dr_time,dr_lat,dr_lon,dr_posErr,dr_u,dr_v,dr_velErr)
        n+=1 
#        print n 
    return drifters
###################################################################################################
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

##################################################################################
def kriging(st,et,sample_step=5,skip=1,num_restarts=10,save=0,output='rbfModel'):
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

   Nd = np.size(xob,1)/skip # number of drifters used
   if (sample_step < 0)|(skip>1): # get samples per time step
      ss = np.abs(sample_step)
      sample_step = 1
      samples = np.arange(0,np.size(tob,0),ss)
      if ss > 1:
         testt = set(np.arange(0,np.size(tob,0))) - set(samples)
         testt = np.array(list(testt))
      else:
         testt = samples
      if skip>1:
         testd = set(np.arange(0,np.size(tob,1))) - set(np.arange(0,np.size(tob,1),skip))
         testd = np.array(list(testd))
      else:
         testd = np.arange(0,np.size(tob,1)) 
      to = np.reshape(tob[samples,::skip],[-1,1])
      yo = np.reshape(yob[samples,::skip],[-1,1])
      xo = np.reshape(xob[samples,::skip],[-1,1])
      uo = np.reshape(yob[samples,::skip],[-1,1])
      vo = np.reshape(xob[samples,::skip],[-1,1])
      tt = np.reshape(tob[testt[:,None],testd],[-1,1])
      yt = np.reshape(yob[testt[:,None],testd],[-1,1])
      xt = np.reshape(xob[testt[:,None],testd],[-1,1])
      ut = np.reshape(yob[testt[:,None],testd],[-1,1])
      vt = np.reshape(xob[testt[:,None],testd],[-1,1])

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


   output_mat = output + '_n'+str(np.size(xo))+'.mat'
   if (ss > 0):
      output_obj_u = output + '_dt'+str(ss*15)+ '_Nd'+str(Nd)+'_u.pkl'
      output_obj_v = output + '_dt'+str(ss*15)+ '_Nd'+str(Nd)+'_v.pkl'
   else:
      output_obj_u = output + '_n'+str(np.size(xo))+'_u.pkl'
      output_obj_v = output + '_n'+str(np.size(xo))+'_v.pkl'
# From here on, always use T,Y,X order
   X = np.concatenate([to,yo,xo],axis=1)
   obs = np.concatenate([vo,uo],axis=1)
   Xt = np.concatenate([tt,yt,xt],axis=1)
   obst = np.concatenate([vt,ut],axis=1)
   rx = 0.1
   sigx = 5
   ry = 0.1
   sigy = 5
   rt = 1.
   sigt = 5
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
   model_v.optimize_restarts(num_restarts=num_restarts)
   print model_v
#   vg,vgVar = model_v.predict(Xrg)
   vHP = model_v.param_array

## Compute U
   ku = k.copy()
   model_u = GPy.models.GPRegression(X,uo,ku)
#   model_u.Gaussian_noise = noise # got from previous experiments
#   model_u.optimize()   
#   print model_u
   model_u.optimize_restarts(num_restarts=num_restarts)
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
      with open(output_obj_u,'wb') as output:
           pickle.dump(model_u,output,-1)
      with open(output_obj_v,'wb') as output:
           pickle.dump(model_v,output,-1)

   print 'End of script, time : ' + str(datetime.now()-startTime)
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
######################################################################################
def getGrid(to,yo,xo,dt=0.5,dx=0.5,xLim=40,yLim=40):
# GRID check size of the final matrix Xrg (reshaped grid)
   if (xo.max()-xo.min())>xLim:
      xmin = xo.mean() - xLim/2
      xmax = xo.mean() + xLim/2
      print 'here x'
   else:
      xmin = xo.min()- dx
      xmax = xo.max() + dx
 
   if (yo.max()-yo.min())>yLim:
      ymin = yo.mean() - yLim/2
      ymax = yo.mean() + yLim/2
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
   return np.concatenate([Tr,Yr,Xr],axis=1)


###################################################################################################
def coKriging(st,et,sample_step=15,nonDiv=1,num_restarts=50,save=0,output='nDivModel'):
   startTime = datetime.now()
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


   if sample_step < 0: # get samples per time step
      ss = sample_step * (-1)
      sample_step = 1
      samples = np.arange(0,np.size(tob,0),ss)
      test = set(np.arange(0,np.size(tob,0))) - set(samples)
      test = np.array(list(test))
      to = np.reshape(tob[samples,:],[-1,1])
      yo = np.reshape(yob[samples,:],[-1,1])
      xo = np.reshape(xob[samples,:],[-1,1])
      uo = np.reshape(yob[samples,:],[-1,1])
      vo = np.reshape(xob[samples,:],[-1,1])
      tt = np.reshape(tob[test,:],[-1,1])
      yt = np.reshape(yob[test,:],[-1,1])
      xt = np.reshape(xob[test,:],[-1,1])
      ut = np.reshape(yob[test,:],[-1,1])
      vt = np.reshape(xob[test,:],[-1,1])

   else:
      ss = 1
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
   obs = np.concatenate([vo,uo],axis=0)
   validPoints = np.where((~np.isnan(xt))&(~np.isnan(yt)))
   tt = tt[validPoints][:,None]
   xt = xt[validPoints][:,None]
   yt = yt[validPoints][:,None]
   ut = ut[validPoints][:,None]
   vt = vt[validPoints][:,None]
   obst = np.concatenate([vt,ut],axis=0)

   output_mat = output + '_n'+str(np.size(xo))+'.mat'
   if ss > 1:
      output_obj = output + '_dt'+str(ss*15)+'.pkl'
   else:
      output_obj = output + '_n'+str(np.size(xo))+'.pkl'

# From here on, always use T,Y,X order
   X = np.concatenate([to,yo,xo],axis=1)
   Xt = np.concatenate([tt,yt,xt],axis=1)

   r = 0.1
   ratio = 0.5
   rt = 1.
   sigt = 5
   noise = 0.0002
#########
# Compute covariances
# Pay attention on the order T,Y,X
   kt = myKernel.Kt(input_dim=1, active_dims=[0], variance=sigt, lengthscale=rt)
   kxy = myKernel.nonDivK(2, [1,2], r)
   k = kt * kxy 
   print k
   model = GPy.models.GPRegression(X,obs,k)
   model.optimize_restarts(num_restarts=num_restarts)
   print model
#   f,fVar = model.predict(Xrg)
#   HP = model.param_array
#   vg = np.reshape(vg,[tg.size,yg.size,-1])
#   ug = np.reshape(ug,[tg.size,yg.size,-1])
#   vgVar = np.reshape(vgVar,[tg.size,yg.size,-1])
#   ugVar = np.reshape(ugVar,[tg.size,yg.size,-1])
   if save==0:
      return X,obs,Xt,obst,model
   else:
      sio.savemat(output_mat,{'Xo':X,'obs':obs,'Xt':Xt,'test_points':obst})
      with open(output_obj,'wb') as output:
           pickle.dump(model,output,-1)
   print 'End of script, time : ' + str(datetime.now()-startTime)
###########################################################################################
def testModel1D(model,Xt,test_points):
    f,fVar = model.predict(Xt)
    return gps.rmse1(f,test_points)
###########################################################################################
def testModel2D(model,Xt,test_points):
    f,fVar = model.predict(Xt)
    return gps.rmse1(f[:f.size/2],test_points[:f.size/2]),rmse1(f[f.size/2:],test_points[f.size/2:])

###########################################################################################
def animateVectors(tg,xg,yg,ug,ugVar,uHP,vg,vgVar,vHP,fname = 'vel_laser.mp4'):

   cbarText = 'Posterior Variance'
   sc = 15
   FS = 30
   figW = 22.
   figH = 12.
   FS = 42

   nframes = np.size(tg)
   X,Y = np.meshgrid(xg,yg)
   U = ug[0,:,:]
   V = vg[0,:,:]
   Var = ugVar[0,:,:]
   Vmin = ugVar.min()
   Vmax = ugVar.max()
# Start Figure
   fig = plt.figure(figsize=(figW,figH))
   plot=fig.add_subplot(111)
   Q = plot.quiver(X, Y, U, V,Var, clim=[Vmin,Vmax], pivot='mid', scale=sc)
   date_text = plot.text(0.02, 1.02, '', transform=plot.transAxes,fontsize=FS,fontweight='bold')

   hpt1 = plot.text(0.4, 1.06, '', transform=plot.transAxes,fontsize=20,fontweight='bold')
   hpt2 = plot.text(0.4, 1.02, '', transform=plot.transAxes,fontsize=20,fontweight='bold')
   hpx1 = plot.text(0.6, 1.06, '', transform=plot.transAxes,fontsize=20,fontweight='bold')
   hpx2 = plot.text(0.6, 1.02, '', transform=plot.transAxes,fontsize=20,fontweight='bold')
   hpy1 = plot.text(0.8, 1.06, '', transform=plot.transAxes,fontsize=20,fontweight='bold')
   hpy2 = plot.text(0.8, 1.02, '', transform=plot.transAxes,fontsize=20,fontweight='bold')
   hpn = plot.text(1., 1.02, '', transform=plot.transAxes,fontsize=20,fontweight='bold')

   hpt1.set_text('$l_t: '+str(np.round(uHP[1],decimals=2))+' h$')
   hpt2.set_text('$\sigma_t: '+str(np.round(uHP[0],decimals=2))+'$')
   hpy1.set_text('$l_y: '+str(np.round(uHP[3],decimals=2))+' km$')
   hpy2.set_text('$\sigma_y: '+str(np.round(uHP[2],decimals=2))+'$')
   hpx1.set_text('$l_x: '+str(np.round(uHP[5],decimals=2))+' km$')
   hpx2.set_text('$\sigma_x: '+str(np.round(uHP[4],decimals=2))+'$')
   hpn.set_text('$noise: '+str(np.round(uHP[6],decimals=2))+'$')

   qk = plot.quiverkey(Q,0.1,0.1,0.5,'0.5m/s')




   plot.tick_params(axis='both',which='major',labelsize=FS)
   plot.set_xlim(xg.min(), xg.max())
   plot.set_ylim(yg.min(), yg.max())
   plot.set_xlabel('Zonal distance (Km)',fontsize=FS)
   plot.set_ylabel('Meridional distance (Km)',fontsize=FS)

   cb = plt.colorbar(Q)
   Vs = np.round(np.arange(Vmin,Vmax,(Vmax-Vmin)/5.),decimals=3)
   cb.set_ticks(Vs)
   cb.ax.tick_params(labelsize = FS)
   cb.set_label(cbarText,fontsize = FS)


   def update_quiver(t, Q, ug, vg,ugVar):
       """updates the horizontal and vertical vector components by a
       fixed increment on each frame
       """
       U = ug[t,:,:]
       V = vg[t,:,:]
       Var = ugVar[t,:,:]
       Q.set_UVC(U,V,Var)
       date_text.set_text('Time: '+ str(tg[t])+' h')
       return Q,
   # you need to set blit=False, or the first set of arrows never gets
   # cleared on subsequent frames
   anim = animation.FuncAnimation(fig, update_quiver, fargs=(Q, ug, vg, ugVar),
                               frames=nframes,interval=5, blit=False)
   
   fps=4
   frame_prefix='_tmp'
   anim.save(fname, fps=fps, codec='mpeg4')
   plt.show()
   return anim

#########################################################################################
def animVecVar(tg,xg,yg,ug,ugVar,uHP,vg,vgVar,vHP):

   sc = 20
   stp = 2
   nframes = np.size(tg)
   X,Y = np.meshgrid(xg,yg)
   U = ug[0,:,:]
   V = vg[0,:,:]
   Uv = ugVar[0,:,:]
   Vv = vgVar[0,:,:]
## gridcolor
   cbarText = 'Posterior Var.'
   cname = 'jet'
   Vmin = np.min([Uv.min(),Vv.min()])  
   Vmax = np.min([Uv.max(),Vv.max()]) 
   N = 30. #number of contour lines
   dL   = (Vmax-Vmin)/N
   Vs = np.arange(Vmin,Vmax+dL,dL) #colorscale
# Figure parameters
   nframes = np.size(tg)
   figW = 22.
   figH = 12.
   FS = 42
# Start Figure
   fig = plt.figure(figsize=(figW,figH))
# start contourf var U
   plot1=fig.add_subplot(231)
   CS1 = plot1.contourf(X,Y,Uv,Vs,cmap=plt.get_cmap(cname))
   plot1.set_xlim(xg.min(), xg.max())
   plot1.set_ylim(yg.min(), yg.max())
   title1 =  plot1.text(0.02, 1.1, '', transform=plot1.transAxes,fontsize=FS,fontweight='bold')
   hpu_text1 = plot1.text(0.02, 0.9, '', transform=plot1.transAxes,fontsize=FS,fontweight='bold')
   hpu_text2 = plot1.text(0.02, 0.8, '', transform=plot1.transAxes,fontsize=FS,fontweight='bold')
   hpu_text3 = plot1.text(0.02, 0.7, '', transform=plot1.transAxes,fontsize=FS,fontweight='bold')
   hpu_text4 = plot1.text(0.02, 0.6, '', transform=plot1.transAxes,fontsize=FS,fontweight='bold')

# start quiver
   plot2=fig.add_subplot(122)
   date_text = plot2.text(0.02, 1.02, '', transform=plot2.transAxes,fontsize=FS,fontweight='bold')
   Q = plot2.quiver(X, Y, U, V, pivot='mid', color='w', scale=sc)
   qk = plot2.quiverkey(Q,0.1,0.1,0.5,'0.5m/s')
   plot2.set_xlim(xg.min(), xg.max())
   plot2.set_ylim(yg.min(), yg.max())

# start contourf var V
   plot3=fig.add_subplot(234)
   CS3 = plot3.contourf(X,Y,Vv,Vs,cmap=plt.get_cmap(cname))
   plot3.set_xlim(xg.min(), xg.max())
   plot3.set_ylim(yg.min(), yg.max())
   title3 =  plot3.text(0.02, 1.1, '', transform=plot3.transAxes,fontsize=FS,fontweight='bold')
   hpv_text1 = plot3.text(0.02, 0.9, '', transform=plot3.transAxes,fontsize=FS,fontweight='bold')
   hpv_text2 = plot3.text(0.02, 0.8, '', transform=plot3.transAxes,fontsize=FS,fontweight='bold')
   hpv_text3 = plot3.text(0.02, 0.7, '', transform=plot3.transAxes,fontsize=FS,fontweight='bold')
   hpv_text4 = plot3.text(0.02, 0.6, '', transform=plot3.transAxes,fontsize=FS,fontweight='bold')

   def anim(j): #, Q, ug, vg):
      global CS1,CS2, Q

      Uv = np.squeeze(ugVar[j,:,:])
      Vv = np.squeeze(vgVar[j,:,:])
      U = ug[j,:,:]
      V = vg[j,:,:]
      CS1 = plot1.contourf(X,Y,Uv,Vs,cmap=plt.get_cmap(cname))
#      cbar1 = plt.colorbar(CS1)
#      cbar1.set_label(cbarText,fontsize = FS)
#      cbar1.ax.tick_params(labelsize = FS)

#      Q = plot.quiver(X,Y,U,V,scale=sc,color='k')    
      Q.set_UVC(U,V)
      date_text.set_text('Time: '+ str(tg[j]))

      CS3 = plot3.contourf(X,Y,Vv,Vs,cmap=plt.get_cmap(cname))
#      cbar3 = plt.colorbar(CS3)
#      cbar1.set_label(cbarText,fontsize = FS)
#      cbar3.ax.tick_params(labelsize = FS)

      print nframes,j
#      return CS,Q
   anim = animation.FuncAnimation(fig, anim,#fargs=(CS,Q),  
                                  frames=nframes, interval=5, blit=False)
   fname = 'velvar_laser.mp4'
   fps=4
   frame_prefix='_tmp'
   anim.save(fname, fps=fps, codec='mpeg4')
   plt.show()
   return anim

################################################################################




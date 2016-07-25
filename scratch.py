import numpy as np
from datetime import datetime


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

def constructField(st,et):

   with open('Filtered_2016_2_7.pkl','rb') as input:
       tr = pickle.load(input)

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
   xob[np.where(np.isnan(lont))]=np.nan
   yob[np.where(np.isnan(lont))]=np.nan

   to = np.repeat(time,np.size(latt,1),axis=1)
   to = np.reshape(to,[-1,1])
   xo = np.reshape(xob,[-1,1])/1000. # in km
   yo = np.reshape(yob,[-1,1])/1000. # in km
   uo = np.reshape(uob,[-1,1])
   vo = np.reshape(vob,[-1,1])

   validPoints = np.where((~np.isnan(xo))&(~np.isnan(yo)))
   uo = uo[validPoints]
   vo = vo[validPoints]
   xo = xo[validPoints]
   yo = yo[validPoints] 
   to = to[validPoints] 
   xo = xo-xo.min() + 2
   yo = yo-yo.min() + 2
   X = np.concatenate([to,xo,yo],axis=1)

   rx = 0.1
   sigx = 5
   ry = 0.1
   sigy = 5
   rt = 1.
   sigy = 5
   noise = 0.0002
#########
# GRID check size of the final matrix Xrg (reshaped grid)
   dt = 0.5
   dx = 0.5
   xg = np.arange(xo.min()-5,xo.max()+5,dx)
   yg = np.arange(yo.min()-5,yo.max()+5,dx)
   tg = np.arange(to.min(),to.max(),dt)
   Tg,Xg,Yg = np.meshgrid(tg,xg,yg) # check if this works
   Tr = np.reshape(Tg,[Tg.size])
   Xr = np.reshape(Xg,[Xg.size])
   Yr = np.reshape(Yg,[Yg.size]) 
   Xrg = np.concatenate(Tr,Xr,Yr)
# Compute covariances

   kt = GPy.kern.RBF(input_dim=1, active_dims=[0], variance=sigt, lengthscale=rt)
   kx = GPy.kern.RBF(input_dim=1, active_dims=[1], variance=sigx, lengthscale=rx)
   ky = GPy.kern.RBF(input_dim=1, active_dims=[2], variance=sigy, lengthscale=ry)
   k = kt * kx * ky 
## Compute U
   model_u = GPy.models.GPRegression(X,uo,k)
#   model_u.Gaussian_noise = noise # got from previous experiments
   model_u.optimize()   
   model_u.optimize_restarts()
   ug,ugVar = model_u.predict(Xrg)
   uHP = model_u.param_array
##
   model_v = GPy.models.GPRegression(X,y2,k)
#   model_v.Gaussian_noise = noise # got from previous experiments
   model_v.optimize()   
   model_u.optimize_restarts()
   vg,vgVar = model_v.predict(Xrg)
   vHP = model_v.param_array

   return tg,xg,yg,ug,ugVar,uHP,vg,vgVar,vHP




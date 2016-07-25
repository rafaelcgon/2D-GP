import numpy as np
from datetime import datetime
from laser_class import drifter,interpolated_tracks,laser
import cPickle as pickle
import os
from dateutil.parser import parse
from datetime import datetime,timedelta
from subprocess import Popen, PIPE
from scipy import interpolate
import GPy
#from joblib import Parallel, delayed  
#import multiprocessing
################################################################
def save_object(obj, filename):
    if np.size(obj) == 1:
       with open(filename, 'wb') as output:
           pickle.dump(obj, output, -1)
    else:
       with open(filename, 'wb') as output:
           for i in range(np.size(obj)):
               pickle.dump(obj[i], output, -1)
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
##################################################################
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
def readFilteredTracks(initdate = datetime(2016,02,07,0,0,0)): # laser_io_methods.py
# Read trajectories from drifters, create objects for each drifter and save all in output
# ID  Date Time   Latitude  Longitude  Pos Error  U (+E-W) V (+N-S) vel Error
# L_0004 2016-01-21 18:45:00.216001  29.03776831107 -87.68717800671    10.9 0.016 0.389 0.033
    
    output = 'Filtered'+'_'+str(initdate.year)+'_'+str(initdate.month)+'_'+str(initdate.day)+'.pkl'
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
    laser_obj = laser(drifters)
    save_object(laser_obj,output)
#    return drifters
################################################################################################
def readDrogueStatus(filename='raw_data/DrogueStatus-Jun11_edited.dat'):
#col1: drifter ID
#col2: drogue status (0=never launched, 1=drogued, 2=undrogued before March 29, 3=unknown)
#col3: date of drogue loss (15-mins iterations since Jan 18).
#col4,5,6,7: date of drogue loss (mon,day,hour,mins).
#col8: launch type (1=LSS, 2=P1, 3=P2, 4=LDA, 5=drogue tests). 

    f = open(filename) 

    drifter_id = []
    drogue_status = []
    drogue_loss_time = []
    drogue_loss_date = []
    launchType = []
    for line in f.readlines():
       if (line[0]!='#'):
          values = [value for value in line.split('  ')]
          if int(values[0])<10:
             id_comp = 'L_000' 
          elif int(values[0])<100:
             id_comp = 'L_00'
          elif int(values[0])<1000:
             id_comp = 'L_0'
          else:
             id_comp = 'L_'
          drifter_id.append(id_comp+str(int(values[0])))
          drogue_status.append(int(values[1]))
          drogue_loss_time.append(values[2])
          mo = np.int(values[3])
          #print mo, values
          da = np.int(values[4])
          hr = np.int(values[5])
          mi = np.int(values[6])
          if mo==0:
             drogue_loss_date.append(datetime(2016,07,01,0,0))
          else:
             drogue_loss_date.append(datetime(2016,mo,da,hr,mi))
          launchType.append(np.int(values[7]))
    return np.array(drifter_id), np.array(drogue_status), np.array(drogue_loss_time), np.array(drogue_loss_date), np.array(launchType)
##########################################################################################
def readLog(date,ship='M8', split = ','):
# read drop positions  
# date: day (and time) when drifters were dropped
#      if date is a vector, it contains the beginnig and end of the period to be analysed  
    if ship == 'M8':
       filename = '/home/rgoncalves/LagOil/LASER/raw_data/M8_drop_20160207.csv'
    elif ship == 'M82':
       filename = '/home/rgoncalves/LagOil/LASER/raw_data/M8_drop_20160211.csv'
    else:
       filename='/home/rgoncalves/LagOil/LASER/raw_data/WS_drifter_drops.csv'
    
    f = open(filename) 
    n = 0
    diff=0
    Matrix=[]
    for line in f.readlines():
#        print line
        values = [value for value in line.split(split)]
        if n>=1:
           if n==1:
              listLen = np.size(values)
              if (ship=='M82'):
                 listLen -= 1 
           if np.size(values)>=listLen:        
           # skip header (GPS wpt,ID,Date,Time,Lon,Lat,)
           # each line is like the following
           # 435,L_0689,2016-02-07,02:06:10,-88.398768,28.820325,1
              if (ship=='M8')or(ship=='M82'): # some drifters from M8 were released together
                 diff = np.size(values) - listLen
              print n,diff,values[diff+2]
              yr = np.int(values[diff+2][0:4])
              mo = np.int(values[diff+2][5:7])
              da = np.int(values[diff+2][8:10])
              hr = np.int(values[diff+3][0:2])
              mi = np.int(values[diff+3][3:5])
              se = np.int(values[diff+3][6:8])
              launchTime = datetime(yr,mo,da,hr,mi,se)
              if launchTime >= date:
                 for i in range(diff+1):                   
#                     drifterNumber = np.int(values[i+1][2:])
                     drifterNumber = values[i+1]
                     print drifterNumber
                     Matrix.append([drifterNumber,launchTime])
        n+=1
    return Matrix  
##########################################################################  
def readTracks(ship='M8'):
# Read trajectories from drifters, create objects for each drifter and save all in output
    if ship =='M82':
       date = datetime(2016,02,10,0,0,0)
    else:
       date = datetime(2016,02,07,0,0,0)
    output = ship+'_'+str(date.year)+'_'+str(date.month)+'_'+str(date.day)+'.pkl'
#first: read time of release of drifters in a log file

    if ship == 'ALL':
       log = readLog(date,'M8')
       log2 = readLog(date,'WS')
       for i in range(np.size(log2,0)):
           log.append(log2[i])
    else:
       log = readLog(date,ship) 
    # sort drifters by name
    log = sorted(log,key=lambda x: x[0])
    # sort drifters by release time
    log = sorted(log,key=lambda x: x[1])
    # number of drifters 
    N_drifters = np.size(log)/2 
    f = open('/home/rgoncalves/LagOil/LASER/raw_data/laser_spot_xml_DRIFTER.dat')
    lines = f.readlines()
    dr=[] # list of drifter objects

    dr_id,dro_stat,dro_loss_time,dro_loss_date,lType = readDrogueStatus()

    for n in range(N_drifters):
        drifter_id = 'L_'+log[n][0][-4:] #just get the number to search the drifter
        print str(n+1)+'/'+str(N_drifters),drifter_id,np.where(dr_id==drifter_id)
        drogueLossDate = dro_loss_date[np.where(dr_id==drifter_id)][0]
        launchType = lType[np.where(dr_id==drifter_id)][0]
        drogueStat = dro_stat[np.where(dr_id==drifter_id)][0]
        
        drifter_data = [s for s in lines if drifter_id in s]
        drifter_id = log[n][0] # get whole drifter id as is in the log
                               # the id letters might be different in the track file
        # get track (time,lat,lon)
        time = [] # date/time when positions are recorder
        time_s = [] # time in seconds (GMT time)
        lat = [] # positions
        lon = [] 
        for line in drifter_data:
            values = [value for value in line.split(' ')]
#'500059883 0-2545404 L_0689 2016-02-03T11:51:50.000Z 1454500310 30.35586 -89.0929 UNLIMITED-TRACK LOW\n'
            yr = np.int(values[3][0:4])
            mo = np.int(values[3][5:7])
            da = np.int(values[3][8:10])
            hr = np.int(values[3][11:13])
            mi = np.int(values[3][14:16])
            se = np.int(values[3][17:19])
            time_aux =datetime(yr,mo,da,hr,mi,se)
            if time_aux >= log[n][1]:
               time.append(time_aux)
               time_s.append(np.int(values[4])) 
               lat.append(np.float(values[5]))
               lon.append(np.float(values[6]))
        # create and save drifter object
        if time != []:
           dr.append(drifter(drifter_id,time,time_s,lat,lon,drogueLossDate,drogueStat,launchType))    
    save_object(dr,output)

########################################################################################## 
def countDataPoints(time,tdata):
    N = np.zeros(np.size(time)-1)
    dt = np.zeros(np.size(time)-1)
    for i in range(np.size(time)-1):
        it = np.where((tdata>=time[i])&(tdata<time[i+1]))
        N[i] = np.size(it)
        dt[i] = np.mean(tdata[it[0][1:]]-tdata[it[0][0:-1]])
    return N,dt

##########################################################################################     
def interp_spline(filename='ALL_2016_2_7.pkl',dt=900,period=10,sp = 0.001):
    # dt is the time step of the interpolated data
    # period is the total period in days of the data to be interpolated
    # sp i the smoothing parameter for the cubic spline interpolation

    output = 'spline_'+filename
    data = read_object('/home/rgoncalves/LagOil/LASER/'+filename)
    # the first drifter released should be data[0]

    period_end = data[0].time[0] + period*86400 # last time step
    time = np.arange(data[0].time[1],period_end,dt)
    time2 = (time-data[0].time[0])/3600.
#    it = np.where((time>=dr.time[0])&(time<=dr.time[-1]))

    N = np.size(data)
    lon = np.zeros((N,np.size(time))) + np.nan
    lat = np.zeros((N,np.size(time))) + np.nan
    u = np.zeros((N,np.size(time))) + np.nan
    v = np.zeros((N,np.size(time))) + np.nan
    drog_stat = np.zeros((N,np.size(time))) # drogue stat
    lastDrogTime = np.zeros(N)  # last time with drogue 
    dLossDate = []  # drogue loss date, in datetime format
    drog_stat0 =  np.zeros(N)  # drogue stat (according to original drogue status file)
    launchType = np.zeros(N) # lauch type
    dr_points = np.zeros((N,np.size(time))) + np.nan # number of data points between interp timesteps
    dr_mdt = np.zeros((N,np.size(time))) + np.nan # mean data frequency between interp timesteps
    dr_id = []
    dr_date0 = []
    dr_time0 = [] 


    for n in range(N):  
        print n
        dr = data[n]
        it = np.where((time>=dr.time[0])&(time<=dr.time[-1]))
        it2 = np.where((dr.time>=time[it][0])&(dr.time<=time[it][-1]))
        it3 = it[0][1:]
#        print it3[1]
        # cubic spline interpolation
        tck_lon = interpolate.splrep(np.array(dr.time)[it2],np.array(dr.lon)[it2],s=sp)
        tck_lat = interpolate.splrep(np.array(dr.time)[it2],np.array(dr.lat)[it2],s=sp)
        lon[n,it] = interpolate.splev(time[it],tck_lon,der=0)
        lat[n,it] = interpolate.splev(time[it],tck_lat,der=0)
        # identify when drogue was lost
        dl = np.where(np.array(dr.date_time)<dr.drogueLoss)
        if np.size(dl) > 0:
           lastDrogTime[n] = dr.time[np.max(dl)]
           print 'last drog time =', lastDrogTime[n]
           print 'last time step =', time[-1]
           drog_stat[n,np.where(time<=lastDrogTime[n])]=1
        else:
           lastDrogTime[n] = -1
           drog_stat[n,:] = 1
        drog_stat0[n] = dr.drogueStat
        dLossDate.append(dr.drogueLoss)
        launchType[n] = dr.launchType
        # Velocity components
        u[n,it] = interpolate.splev(time[it],tck_lon,der=1)*111000*np.cos(lat[n,it]*np.pi/180.)
        v[n,it] = interpolate.splev(time[it],tck_lat,der=1)*111000
        M1,M2 = countDataPoints(time[it],np.array(dr.time))
        dr_points[n,it3] = M1
        dr_points[n,it3[0]-1] = M1[0]
        dr_mdt[n,it3] = M2
        dr_mdt[n,it3[0]-1] = M2[0]
        dr_id.append(dr.id)
        dr_date0.append(dr.date_time[0])
        dr_time0.append(dr.time[0])
    itracks = interpolated_tracks(dr_id,time2,lon,lat,u,v,dr_date0,dr_time0,dr_points,dr_mdt, 
                drog_stat0,drog_stat,lastDrogTime,dLossDate,launchType)
    save_object(itracks,output)
#    return itracks
##########################################################################################     
def interp_kriging(filename='ALL_2016_2_7.pkl',dt=900,period=10,checkError=0):
    # dt is the time step of the interpolated data
    # period is the total period in days of the data to be interpolated
    # sp i the smoothing parameter for the cubic spline interpolation

    output = 'kriging2_'+filename
    data = read_object('/home/rgoncalves/LagOil/LASER/'+filename)
    if checkError ==1:
       with open('spline_'+filename,'rb') as input:
            spl = pickle.load(input)
       
    # the first drifter released should be data[0]

    period_end = data[0].time[0] + period*86400 # last time step
    time = np.arange(data[0].time[1],period_end,dt)
    time2 = (time-data[0].time[0])/3600. #time in hours 
    time2 = np.reshape(time2,[time2.size,1])
    timem= (time2[1:]+time2[0:-1])/2.
#    it = np.where((time>=dr.time[0])&(time<=dr.time[-1]))

    N = np.size(data)
    lon = np.zeros((N,np.size(time))) + np.nan
    lat = np.zeros((N,np.size(time))) + np.nan
    varLon = np.zeros((N,np.size(time))) + np.nan
    varLat = np.zeros((N,np.size(time))) + np.nan
    lenLon = np.zeros(N)
    lenLat = np.zeros(N)
    varianceLon = np.zeros(N)
    varianceLat = np.zeros(N)
    noiseLon = np.zeros(N)
    noiseLat = np.zeros(N)
#    variance = np.zeros((N,np.size(time))) + np.nan
    lon2 = np.zeros((N,np.size(time)-1)) + np.nan
    lat2 = np.zeros((N,np.size(time)-1)) + np.nan
    u = np.zeros((N,np.size(time)-1)) + np.nan
    v = np.zeros((N,np.size(time)-1)) + np.nan
# drogue stat
    drog_stat = np.zeros((N,np.size(time))) 
# last time with drogue
    lastDrogTime = np.zeros(N)   
# drogue loss date, in datetime format
    dLossDate = []  
# drogue stat (according to original drogue status file)
    drog_stat0 =  np.zeros(N)  
# lauch type
    launchType = np.zeros(N) 
# number of data points between interp timesteps
    dr_points = np.zeros((N,np.size(time)-1)) + np.nan 
# mean data frequency between interp timesteps
    dr_mdt = np.zeros((N,np.size(time)-1)) + np.nan 
    dr_id = []
    dr_date0 = []
    dr_time0 = [] 
# Loop over every drifter
    kernel = GPy.kern.RBF(input_dim=1, variance=1159.68, lengthscale=4.5)
    for n in range(N):  
        print n
        print kernel
        dr = data[n]
        it = np.where((time>=dr.time[0])&(time<=dr.time[-1]))
        it2 = np.where((dr.time>=time[it][0])&(dr.time<=time[it][-1]))
        it3 = it[0][1:]
# gaussian process interpolation
        X = np.array(dr.time)[it2][:,None]
        Y1 = np.array(dr.lon)[it2][:,None]
        Y2 = np.array(dr.lat)[it2][:,None]
        if (it2[0][0] > 0):
          X = np.concatenate([np.array([dr.time[it2[0][0]-1]])[:,None],X],axis=0)
          Y1 = np.concatenate([np.array([dr.lon[it2[0][0]-1]])[:,None],Y1],axis=0)
          Y2 = np.concatenate([np.array([dr.lat[it2[0][0]-1]])[:,None],Y2],axis=0)  
        if (dr.time[-1] > X[-1]):
          X = np.concatenate([X,np.array([dr.time[it2[0][-1]+1]])[:,None]],axis=0)
          Y1 = np.concatenate([Y1,np.array([dr.lon[it2[0][-1]+1]])[:,None]],axis=0)
          Y2 = np.concatenate([Y2,np.array([dr.lat[it2[0][-1]+1]])[:,None]],axis=0)         

        X = (X-data[0].time[0])/3600.
#        Y = np.array([Y1,Y2]).T
        model_lon = GPy.models.GPRegression(X,Y1,kernel)
        model_lon.Gaussian_noise = 1.75598244486e-07 # got from previous experiments
        variable1,variance1 = model_lon.predict(time2[it])

        model_lat = GPy.models.GPRegression(X,Y2,kernel)
        model_lat.Gaussian_noise = 1.75598244486e-07 # got from previous experiments
        variable2,variance2= model_lat.predict(time2[it])   

        if checkError==0:
           model_lon.optimize(messages=True)
           variable1,vari1 = model_lon.predict(time2[it]) 
           model_lat.optimize(messages=True)
           variable2,vari2 = model_lat.predict(time2[it])       
        elif checkError==1:
           error_lo = np.square(variable1-spl.lon[n,it])
           error_la = np.square(lat[variable2,:]-spl.lat[n,it])
           rmse_lon = np.sqrt(np.nanmean(error_lo))
           rmse_lat = np.sqrt(np.nanmean(error_la))
           if (rmse_lon>1e-3):           
              model_lon.optimize(messages=True)
              variable1,vari1 = model_lon.predict(time2[it]) 
              print model_lon
           if (rmse_lat>1e-3):
              model_lat.optimize(messages=True)
              variable2,vari2 = model_lat.predict(time2[it]) 
        lon[n,it] = variable1.T 
        varLon[n,it] = variance1.T
        lenLon[n] = model_lon.rbf.lengthscale[0]
        varianceLon[n] = model_lon.rbf.variance[0]
        noiseLon[n] = model_lon.Gaussian_noise[0]

        lat[n,it] = variable2.T 
        varLat[n,it] = variance2.T 
        lenLat[n] = model_lat.rbf.lengthscale[0]
        varianceLat[n] = model_lat.rbf.variance[0]
        noiseLat[n] = model_lat.Gaussian_noise[0]

#           kernel = GPy.kern.RBF(input_dim=1, variance=1159.68, lengthscale=4.5)
#        variance[n,it] = vari.T
# variance is sigma**2!
# variance * np.exp(-np.square(x-xo)/(2*np.square(lengthscale))))
#
# for spline interpolation:
#        tck_lon = interpolate.splrep(np.array(dr.time)[it2],np.array(dr.lon)[it2],s=sp)
#        tck_lat = interpolate.splrep(np.array(dr.time)[it2],np.array(dr.lat)[it2],s=sp)
#        lon[n,it] = interpolate.splev(time[it],tck_lon,der=0)
#        lat[n,it] = interpolate.splev(time[it],tck_lat,der=0)

        # identify when drogue was lost
        dl = np.where(np.array(dr.date_time)<dr.drogueLoss)
        if np.size(dl) > 0:
           lastDrogTime[n] = dr.time[np.max(dl)]
           print 'last drog time =', lastDrogTime[n]
           print 'last time step =', time[-1]
           drog_stat[n,np.where(time<=lastDrogTime[n])]=1
        else:
           lastDrogTime[n] = -1
           drog_stat[n,:] = 1
        drog_stat0[n] = dr.drogueStat
        dLossDate.append(dr.drogueLoss)
        launchType[n] = dr.launchType
# Velocity components NOT READY
        timeDiff = np.diff(time[it],n=1)
        latm = (np.squeeze(lat[n,it])[1:]+np.squeeze(lat[n,it])[:-1])/2.
        lonm = (np.squeeze(lon[n,it])[1:]+np.squeeze(lon[n,it])[:-1])/2.
        u[n,it3-1] = np.diff(lon[n,it], n=1, axis=-1)*111000*np.cos(latm*np.pi/180.)/timeDiff
        v[n,it3-1] = np.diff(lat[n,it], n=1, axis=-1)*111000/timeDiff
#        lat2[n,it3-1] = latm
#        lon2[n,it3-1] = lonm
#        u[n,it] = interpolate.splev(time[it],tck_lon,der=1)*111000*np.cos(lat[n,it]*np.pi/180.)
#        v[n,it] = interpolate.splev(time[it],tck_lat,der=1)*111000
        
        M1,M2 = countDataPoints(time[it],np.array(dr.time))
        dr_points[n,it3-1] = M1
#        dr_points[n,it3[0]-1] = M1[0]
        dr_mdt[n,it3-1] = M2
#        dr_mdt[n,it3[0]-1] = M2[0]
        dr_id.append(dr.id)
        dr_date0.append(dr.date_time[0])
        dr_time0.append(dr.time[0])
    itracks = interpolated_tracks(dr_id,time2,lon,lat,u,v,dr_date0,dr_time0,dr_points,dr_mdt, 
                drog_stat0,drog_stat,lastDrogTime,dLossDate,launchType,varLon,varLat, 
                lenLon,lenLat,varianceLon,varianceLat,noiseLon,noiseLat)
    save_object(itracks,output)
#    return lenLon,lenLat,varianceLon,varianceLat,noiseLon,noiseLat



##########################################################################################     
def interp_kriging2(filename='ALL_2016_2_7.pkl',dt=900,period=10):
    # dt is the time step of the interpolated data
    # period is the total period in days of the data to be interpolated
    # sp i the smoothing parameter for the cubic spline interpolation
    output = 'kriging2_'+filename
    data = read_object('/home/rgoncalves/LagOil/LASER/'+filename)

    # the first drifter released should be data[0]

    period_end = data[0].time[0] + period*86400 # last time step
    time = np.arange(data[0].time[1],period_end,dt)
#    time2 = (time-data[0].time[0])/3600. #time in hours 
#    time2 = np.reshape(time2,[time2.size,1])
#    timem= (time2[1:]+time2[0:-1])/2.
#    it = np.where((time>=dr.time[0])&(time<=dr.time[-1]))

    N = np.size(data)
    lon = np.zeros((N,np.size(time))) + np.nan
    lat = np.zeros((N,np.size(time))) + np.nan
    varLon = np.zeros((N,np.size(time))) + np.nan
    varLat = np.zeros((N,np.size(time))) + np.nan
    lenLon = np.zeros(N)
    lenLat = np.zeros(N)
    varianceLon = np.zeros(N)
    varianceLat = np.zeros(N)
    noiseLon = np.zeros(N)
    noiseLat = np.zeros(N)
#    variance = np.zeros((N,np.size(time))) + np.nan
    lon2 = np.zeros((N,np.size(time)-1)) + np.nan
    lat2 = np.zeros((N,np.size(time)-1)) + np.nan
    u = np.zeros((N,np.size(time)-1)) + np.nan
    v = np.zeros((N,np.size(time)-1)) + np.nan
# drogue stat
    drog_stat = np.zeros((N,np.size(time))) 
# last time with drogue
    lastDrogTime = np.zeros(N)   
# drogue loss date, in datetime format
    dLossDate = []  
# drogue stat (according to original drogue status file)
    drog_stat0 =  np.zeros(N)  
# lauch type
    launchType = np.zeros(N) 
# number of data points between interp timesteps
    dr_points = np.zeros((N,np.size(time)-1)) + np.nan 
# mean data frequency between interp timesteps
    dr_mdt = np.zeros((N,np.size(time)-1)) + np.nan 
    dr_id = []
    dr_date0 = []
    dr_time0 = []

 
    num_cores = multiprocessing.cpu_count()
    dr_id, time2, lon, lat, u, v, dr_date0, dr_time0, dr_points, dr_mdt, drog_stat0, drog_stat, lastDrogTime, dLossDate, launchType, varLon, varLat, lenLon, lenLat, varianceLon, varianceLat, noiseLon, noiseLat = Parallel(n_jobs=num_cores)(delayed(kriging)(data[i],time) for i in range(N))

    itracks = interpolated_tracks(dr_id,time2,lon,lat,u,v,dr_date0,dr_time0,dr_points,dr_mdt, 
                drog_stat0,drog_stat,lastDrogTime,dLossDate,launchType,varLon,varLat, 
                lenLon,lenLat,varianceLon,varianceLat,noiseLon,noiseLat)
    save_object(itracks,output)
#    return lenLon,lenLat,varianceLon,varianceLat,noiseLon,noiseLat

#############################################################################################
def kriging(dr,time):
    # initialize arrays
    time2 = time - time[0]
    lon = np.zeros(np.size(time)) + np.nan
    lat = np.zeros(np.size(time)) + np.nan
    u = np.zeros(np.size(time)-1) + np.nan
    v = np.zeros(np.size(time)-1) + np.nan
    drog_stat = np.zeros(np.size(time))
 
    it = np.where((time>=dr.time[0])&(time<=dr.time[-1]))
    it2 = np.where((dr.time>=time[it][0])&(dr.time<=time[it][-1]))
    it3 = it[0][1:]
# gaussian process interpolation
    kernel = GPy.kern.RBF(input_dim=1, variance=1159.68, lengthscale=4.5)
    X = np.array(dr.time)[it2][:,None]
    Y1 = np.array(dr.lon)[it2][:,None]
    Y2 = np.array(dr.lat)[it2][:,None]
    if (it2[0][0] > 0):
       X = np.concatenate([np.array([dr.time[it2[0][0]-1]])[:,None],X],axis=0)
       Y1 = np.concatenate([np.array([dr.lon[it2[0][0]-1]])[:,None],Y1],axis=0)
       Y2 = np.concatenate([np.array([dr.lat[it2[0][0]-1]])[:,None],Y2],axis=0)  
    if (dr.time[-1] > X[-1]):
       X = np.concatenate([X,np.array([dr.time[it2[0][-1]+1]])[:,None]],axis=0)
       Y1 = np.concatenate([Y1,np.array([dr.lon[it2[0][-1]+1]])[:,None]],axis=0)
       Y2 = np.concatenate([Y2,np.array([dr.lat[it2[0][-1]+1]])[:,None]],axis=0)         
    X = (X - time[0])/3600.

    model_lon = GPy.models.GPRegression(X,Y1,kernel)
    model_lon.Gaussian_noise = 1.75598244486e-07 # got from previous experiments
    model_lon.optimize(messages=True)
    lenLon = model_lon.rbf.lengthscale[0]
    varianceLon = model_lon.rbf.variance[0]
    noiseLon = model_lon.Gaussian_noise[0]
    variable,variance = model_lon.predict(time2[it])
    lon[it] = variable 
    varLon[it] = variance
        
    model_lat = GPy.models.GPRegression(X,Y2,kernel)
    model_lat.Gaussian_noise = 1.75598244486e-07 # got from previous experiments
    model_lat.optimize(messages=True)
    lenLat = model_lat.rbf.lengthscale[0]
    varianceLat = model_lat.rbf.variance[0]
    noiseLat = model_lat.Gaussian_noise[0]
    variable,variance= model_lat.predict(time2[it])   
    lat[it] = variable 
    varLat[it] = variance.T 
    # identify when drogue was lost
    dl = np.where(np.array(dr.date_time)<dr.drogueLoss)
    if np.size(dl) > 0:
       lastDrogTime = dr.time[np.max(dl)]
       print 'last drog time =', lastDrogTime
       print 'last time step =', time[-1]
       drog_stat[np.where(time<=lastDrogTime)]=1
    else:
       lastDrogTime = -1
       drog_stat[:] = 1
# Velocity components NOT READY
    timeDiff = np.diff(time[it],n=1)
    latm = (np.squeeze(lat[it])[1:]+np.squeeze(lat[it])[:-1])/2.
    lonm = (np.squeeze(lon[it])[1:]+np.squeeze(lon[it])[:-1])/2.
    u[it3-1] = np.diff(lon[it], n=1, axis=-1)*111000*np.cos(latm*np.pi/180.)/timeDiff
    v[it3-1] = np.diff(lat[it], n=1, axis=-1)*111000/timeDiff
        
    M1,M2 = countDataPoints(time[it],np.array(dr.time))
    dr_points[it3-1] = M1
#        dr_points[n,it3[0]-1] = M1[0]
    dr_mdt[n,it3-1] = M2
#        dr_mdt[n,it3[0]-1] = M2[0]

    return dr.id, time2, lon, lat, u, v, dr.date_time[0], dr.time[0], dr_points, dr_mdt, dr.drogStat, drog_stat, lastDrogTime, dr.drogueLoss, dr.launchType, varLon, varLat, lenLon, lenLat, varianceLon, varianceLat, noiseLon, noiseLat

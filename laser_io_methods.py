import numpy as np
from datetime import datetime
from laser_class import drifter,interpolated_tracks
import cPickle as pickle
import os
from dateutil.parser import parse
from datetime import datetime,timedelta
from subprocess import Popen, PIPE
from scipy import interpolate
import GPy
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

    output = 'interp_'+filename
    data = read_object(filename)
    # the first drifter released should be data[0]

    period_end = data[0].time[0] + period*86400 # last time step
    time = np.arange(data[0].time[0],period_end,dt)

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
    itracks = interpolated_tracks(dr_id,time,lon,lat,u,v,dr_date0,dr_time0,dr_points,dr_mdt, 
                drog_stat0,drog_stat,lastDrogTime,dLossDate,launchType)
    save_object(itracks,output)
#    return itracks
##########################################################################################     
def interp_kriging(filename='ALL_2016_2_7.pkl',dt=900,period=10,sp = 0.001):
    # dt is the time step of the interpolated data
    # period is the total period in days of the data to be interpolated
    # sp i the smoothing parameter for the cubic spline interpolation

    output = 'kriging_'+filename
    data = read_object(filename)
    # the first drifter released should be data[0]

    period_end = data[0].time[0] + period*86400 # last time step
    time = np.arange(data[0].time[0],period_end,dt)

#    it = np.where((time>=dr.time[0])&(time<=dr.time[-1]))

    N = np.size(data)
    lon = np.zeros((N,np.size(time))) + np.nan
    lat = np.zeros((N,np.size(time))) + np.nan
    varlon = np.zeros((N,np.size(time))) + np.nan
    varlat = np.zeros((N,np.size(time))) + np.nan
    u = np.zeros((N,np.size(time))) + np.nan
    v = np.zeros((N,np.size(time))) + np.nan
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
    dr_points = np.zeros((N,np.size(time))) + np.nan 
# mean data frequency between interp timesteps
    dr_mdt = np.zeros((N,np.size(time))) + np.nan 
    dr_id = []
    dr_date0 = []
    dr_time0 = [] 
# Loop over every drifter
    for n in range(N):  
        print n
        dr = data[n]
        it = np.where((time>=dr.time[0])&(time<=dr.time[-1]))
        it2 = np.where((dr.time>=time[it][0])&(dr.time<=time[it][-1]))
        it3 = it[0][1:]
# gaussian process interpolation
        kernel = GPy.kern.RBF(input_dim=1, variance=100., lengthscale=2)
        X = np.array(dr.time)[it2]
        X = np.reshape(X,[X.size,1])
        Y1 = np.array(dr.lon)[it2]
        Y2 = np.array(dr.lat)[it2]
        Y = np.array([Y1,Y2]).T
        model = GPy.models.GPRegression(X,Y,kernel)
        model.optimize()
        variables,variance[n,it] = model.predict(time[it])
#        lon[n,it] = interpolate.splev(time[it],tck_lon,der=0)
#        lat[n,it] = interpolate.splev(time[it],tck_lat,der=0)
        lon[n,it] = variables[:,0]
        lat[n,it] = variables[:,1]   
# variance is sigma**2!
# variance * np.exp(-np.square(x-xo)/(2*np.square(lengthscale))))

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
#        u[n,it] = interpolate.splev(time[it],tck_lon,der=1)*111000*np.cos(lat[n,it]*np.pi/180.)
#        v[n,it] = interpolate.splev(time[it],tck_lat,der=1)*111000
        M1,M2 = countDataPoints(time[it],np.array(dr.time))
        dr_points[n,it3] = M1
        dr_points[n,it3[0]-1] = M1[0]
        dr_mdt[n,it3] = M2
        dr_mdt[n,it3[0]-1] = M2[0]
        dr_id.append(dr.id)
        dr_date0.append(dr.date_time[0])
        dr_time0.append(dr.time[0])
    itracks = interpolated_tracks(dr_id,time,lon,lat,u,v,dr_date0,dr_time0,dr_points,dr_mdt, 
                drog_stat0,drog_stat,lastDrogTime,dLossDate,launchType)
    save_object(itracks,output)
#    return itracks






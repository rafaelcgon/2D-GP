from datetime import datetime
import numpy as np

class drifter(object):
    def __init__(self,drifter_id,date_time,time,lat,lon,droLossDate,droStat,launchType):
        self.id = drifter_id
        self.date_time = date_time
        self.time = time
        self.lat = lat
        self.lon = lon
        self.drogueLoss = droLossDate
        self.drogueStat = droStat
        self.launchType = launchType
#        self.timeStamp = drifter_timeStamp
##########################################################
class interpolated_tracks(object):
    def __init__(self,drifter_id,time,lon,lat,u,v,date0,time0,nsamples,mdt,drog_stat0,
                   drog_stat,lDrogueTime,dLossDate,launchType,varLon=0,varLat=0,
                   lenLon=0,lenLat=0,varianceLon=0,varianceLat=0,noiseLon=0,noiseLat=0):
        self.id = drifter_id      # drifter ID
        self.time = time          # interpolated time, in seconds
        self.date0 = date0        # date and time of first data collected
        self.time0 = time0        # time in seconds of first data collected
        self.lat = lat            # latitude [drifter,time]
        self.lon = lon            # longitude [drifter,time]
        if (np.size(varLon)>1): # uncertainty of pos estimate (posterior variance)
           self.pos_varLon = varLon # use only for kriging
           self.pos_varLat = varLat
       # optimized hyperparameters for each drifter
           self.lenLon = lenLon
           self.lenLat = lenLat
           self.varianceLon = varianceLon
           self.varianceLat = varianceLat
           self.noiseLon = noiseLon
           self.noiseLat = noiseLat
        self.u = u                # zonal velocity 
        self.v = v                # meridional velocity
        self.n_samples = nsamples # number of samples before timestep n after timestep n-1
        self.data_freq = mdt      # average data frequency btwn 2 time steps
        self.drogueStat0 = drog_stat0 # status according to file DrogueStatus-Jun11.dat
        self.drogueStat = drog_stat # time series of the drogue status(0=undrogued, 1=drogued)
        self.lastDrogueTime = lDrogueTime # last time with drogue (==-1 if didn't lose )
        self.lossDate = dLossDate # date when drogue was lost
        self.launchType = launchType # (1=LSS, 2=P1, 3=P2, 4=LDA, 5=drogue tests)
##########################################################
class triangles(object):
    def __init__(self,time,Area,Sides,Nodes_id,Lon_n,Lat_n,
                 Lon_c,Lat_c,U_n,V_n,Div,Curl,Stretch,Shear):
        self.time = np.array([time])
        self.Area = np.array([Area])                # triangle area
        self.Sides = np.reshape(Sides,[3,1])       # sides of the triangle
        self.Nodes_id = np.reshape(Nodes_id,[3,1]) # id of the 3 drifters
        self.Lon_n = np.reshape(Lon_n,[3,1])       # longitude of the 3 drifters
        self.Lat_n = np.reshape(Lat_n,[3,1])       # latitude of the 3 drifters
        self.Un = np.reshape(U_n,[3,1])            # zonal velocity of the 3 drifters
        self.Vn = np.reshape(V_n,[3,1])            # meridional velocity of the 3 drifters
        self.Lon_c = np.array([Lon_c])     # longitude of the centroid
        self.Lat_c = np.array([Lat_c])     # latitude of the centroid
        self.Div = np.array([Div])
        self.Curl = np.array([Curl])
        self.Stretch = np.array([Stretch])
        self.Shear = np.array([Shear])
###
    def addData(self,time,Area,Sides,Lon_n,Lat_n,
                 Lon_c,Lat_c,U_n,V_n,Div,Curl,Stretch,Shear):
        self.time = np.append(self.time,time)
        self.Area = np.append(self.Area,Area)  
        self.Sides = np.append(self.Sides,Sides,axis=1) 
#        self.Nodes_id = np.append(self.Nodes_id,np.reshape(Nodes_id,[1,-1]),axis=0)
        self.Lon_n = np.append(self.Lon_n,Lon_n,axis=1)
        self.Lat_n = np.append(self.Lat_n,Lat_n,axis=1)
        self.Un = np.append(self.Un,U_n,axis=1)
        self.Vn = np.append(self.Vn,V_n,axis=1)
        self.Lon_c = np.append(self.Lon_c,Lon_c)  
        self.Lat_c = np.append(self.Lat_c,Lat_c)  
        self.Div = np.append(self.Div,Div)
        self.Curl = np.append(self.Curl,Curl)
        self.Stretch = np.append(self.Stretch,Stretch)
        self.Shear = np.append(self.Shear,Shear)
##########################################################
class ASCAT(object):
    def __init__(self,time,lat,lon,vwnd,uwnd,wnd_speed,vwnd_rms,uwnd_rms,
                     wnd_speed_rms,samp_length):
        self.time = time
        self.lat = lat
        self.lon = lon
        self.vwnd = vwnd
        self.uwnd = uwnd
        self.wnd_speed = wnd_speed
        self.vwnd_rms = vwnd_rms
        self.uwnd_rms = uwnd_rms
        self.wnd_speed_rms = wnd_speed_rms
        self.samp_length = samp_length




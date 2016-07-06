import numpy as np
import scipy as sc
from scipy.interpolate import Rbf
from matplotlib import pylab as pl
import matplotlib.pyplot as plt               #
from matplotlib import rc,rcParams
import matplotlib as mpl 
from GP_scripts import *
import sys
import cPickle as pickle
sys.path.append('/home/rgoncalves/LagOil/LASER')
#import laser_post_process as lpp
from geopy.distance import vincenty,GreatCircleDistance
import mpl_toolkits.basemap.pyproj as pyproj

def laser(ts=20,nsteps=8,l_df=5,l_cf=5,rate=0.5,noise = 0.0025,nsamples=1):
#
#  Load Observations
   with open('interp_ALL_2016_2_7.pkl','rb') as input:
        tr = pickle.load(input)
    
   #last = 343 # Use this in case you want the initial time to be after the last drifter is launched (for dt = 1800 s)
   # start time index
   st = 0 #np.where(np.abs(tr.time-tr.time0[last])==np.min(np.abs(tr.time-tr.time0[last])))[0][0]+1
   # end time index
   # for dt = 1800 s (1/2 hour), there are 48 samples per day 
   # for dt = 900 s (1/4 hour), there are 96 samples per day
   # compute fsle for 5 days et =  start + 48 * 5
   et = st + (96 * 7) + 1
   # get data for the period to be considered
   time = tr.time[st:et] - tr.time[st]    
   latt = tr.lat[:,st:et]
   lont = tr.lon[:,st:et]
   drogue = tr.drogueStat[:,st:et]
   uob = tr.u[:,st:et]
   vob = tr.v[:,st:et]
   uob[np.where(drogue==0)]=np.nan 
   vob[np.where(drogue==0)]=np.nan 
   # origin of cartesian coord.
   lat0 = 28.8
   lon0 = -88.55
# Put observations in cartesian coord.
#   xob = np.zeros((np.size(lont,0),np.size(lont,1)))
#   yob = np.zeros((np.size(lont,0),np.size(lont,1)))
#   for t in range(np.size(lont,1)):
#      latt2 = np.squeeze(latt[:,t])
#      lont2 = np.squeeze(lont[:,t])
#      for i in range(np.size(lont,0)):
#         hip = GreatCircleDistance((latt2[i],lont2[i]),(lat0,lon0)).km
#         xob[i,t] = GreatCircleDistance((latt2[i],lont2[i]),(latt2[i],lon0)).km
#         yob[i,t] = GreatCircleDistance((latt2[i],lont2[i]),(lat0,lont2[i])).km

   NAD83=pyproj.Proj("+init=EPSG:3452") #Louisiana South (ftUS)
   xob,yob=NAD83(lont,latt)
   xob[np.where(np.isnan(lont))]=np.nan
   yob[np.where(np.isnan(lont))]=np.nan

   xo = np.reshape(xob[:,ts:ts+nsteps],[-1])/1000. # in km
   yo = np.reshape(yob[:,ts:ts+nsteps],[-1])/1000. # in km
   uo = np.reshape(uob[:,ts:ts+nsteps],[-1])
   vo = np.reshape(vob[:,ts:ts+nsteps],[-1])
   uo[np.where(np.abs(uo)>2)]=np.nan 
   vo[np.where(np.abs(vo)>2)]=np.nan

   invPoints = np.where((~np.isnan(uo))&(~np.isnan(vo)))
   uo = uo[invPoints]
   vo = vo[invPoints]
   xo = xo[invPoints]#np.where((state==1)&(~np.isnan(xo)))]
   yo = yo[invPoints] #np.where((state==1)&(~np.isnan(yo)))]
   xo = xo-xo.min() + 2
   yo = yo-yo.min() + 2
   
   if nsamples>0:
      samples = np.arange(0,xo.size,3) #np.random.randint(0,xo.size,nsamples)
      test = set(np.arange(xo.size)) - set(samples)
      test = np.array(list(test))
      xt = xo[test] # use to compute error
      yt = yo[test]
      ut = uo[test]
      vt = vo[test]
      xo = xo[samples]
      yo = yo[samples]
      uo = uo[samples]
      vo = vo[samples]
   else:
      xt = np.array([0])
      yt = np.array([0])
      ut = np.array([0])
      vt = np.array([0])
      
   obs = np.concatenate([uo,vo])
   obs = np.reshape(obs,[obs.size,1])
#########
# GRID
   dx = 0.5
   x = np.arange(np.min([xo.min(),xt.min()])-5,
                 np.max([xo.max(),xt.max()])+5,dx)
   y = np.arange(np.min([yo.min(),yt.min()])-5,
                 np.max([yo.max(),yt.max()])+5,dx)
   X,Y = np.meshgrid(x,y)
   Xs = np.reshape(X,[X.size])
   Ys = np.reshape(Y,[Y.size])   
# Compute covariances
   # obs points
   print 'compute K'
   K = rate*compute_K(xo,yo,l_df,1) + (1-rate)*compute_K(xo,yo,l_cf,2) 
   Ko = np.identity(np.size(K,0))*noise # obs noise
   K = K + Ko 
   # invert K
   print 'invert K'
   Ki = np.linalg.inv(K)
   # test points
   print 'compute Ks'

   Ks = rate*compute_Ks(xo,yo,Xs,Ys,l_df,1) + (1-rate)*compute_Ks(xo,yo,Xs,Ys,l_cf,2)

   # verification points
   print 'compute Kst'
   Kst = rate*compute_Ks(xo,yo,xt,yt,l_df,1) + (1-rate)*compute_Ks(xo,yo,xt,yt,l_cf,2)
   print 'compute Kss'
   Kss = rate*compute_K(Xs,Ys,l_df,1) + (1-rate)*compute_K(Xs,Ys,l_cf,2)
   Cov = Kss - np.dot(Ks,np.dot(Ki,Ks.T))
   uvar =  np.reshape(np.diag(Cov[:X.size,:X.size]),[y.size,-1])
   vvar =  np.reshape(np.diag(Cov[X.size:,X.size:]),[y.size,-1])

# get mean 
   f = getMean(Ks,Ki,obs)
   uf = np.reshape(f[:f.size/2],[y.size,-1])
   vf = np.reshape(f[f.size/2:],[y.size,-1])

   ft = getMean(Kst,Ki,obs)
   uft = ft[:ft.size/2]
   vft = ft[ft.size/2:]

   return x,y,uf,vf,xo,yo,uo,vo,uvar,vvar,xt,yt,ut,vt,uft,vft 

########################################################################################
def simLaser(ts=0,l_df = 2,l_cf = 2,rate=0.5,noise = 0.05):
  
   dx = 0.5
   x = np.arange(0,25+dx,dx)
   y = np.arange(0,25+dx,dx)
   X,Y = np.meshgrid(x,y)
   Xs = np.reshape(X,[X.size])
   Ys = np.reshape(Y,[Y.size])   

   with open('simulTracks.pkl','rb') as input:
        tr = pickle.load(input)

   lat0 = 28.69
   lon0 = -88.28
# Put observations in a 20 km by 20 km dimension
   xob = np.zeros((np.size(tr.lon,0),np.size(tr.lon,1)))
   yob = np.zeros((np.size(tr.lon,0),np.size(tr.lon,1)))
   for t in range(np.size(tr.lon,1)):
      latt = tr.lat[:,t]
      lont = tr.lon[:,t]
      for i in range(np.size(tr.lon,0)):
         xob[i,t] = GreatCircleDistance((latt[i],lont[i]),(latt[i],lon0)).km
         yob[i,t] = GreatCircleDistance((latt[i],lont[i]),(lat0,lont[i])).km
#      xo[i] = vincenty((latt[i],lont[i]),(latt[i],lon0)).km
#      yo[i] = vincenty((latt[i],lont[i]),(lat0,lont[i])).km
   xo = xob[:,ts]
   yo = yob[:,ts]
   uo = tr.u[:,ts]
   vo = tr.v[:,ts]
   obs = np.concatenate([uo,vo])
   obs = np.reshape(obs,[obs.size,1])

   K = rate*compute_K(xo,yo,l_df,1) + (1-rate)*compute_K(xo,yo,l_cf,2) 
   Ko = np.identity(np.size(K,0))*noise
   K = K + Ko 
   Ki = np.linalg.inv(K)
   Ks = rate*compute_Ks(xo,yo,Xs,Ys,l_df,1) + (1-rate)*rate*compute_Ks(xo,yo,Xs,Ys,l_cf,2)
   
   f = getMean(Ks,Ki,obs)
   uf = np.reshape(f[:f.size/2],[y.size,-1])
   vf = np.reshape(f[f.size/2:],[y.size,-1])

   return X,Y,uf,vf,xob,yob,tr.u,tr.v

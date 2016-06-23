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
import laser_post_process as lpp
from geopy.distance import vincenty,GreatCircleDistance

def laser(timestep,sigma = 10):
#
#  Load Observations
   with open('/home/rgoncalves/LagOil/LASER/interp_ALL_2016_2_7.pkl','rb') as input:
        tr = pickle.load(input)

   last = 343 
   # start time index
   st = 0 #np.where(np.abs(tr.time-tr.time0[last])==np.min(np.abs(tr.time-tr.time0[last])))[0][0]+1
   # end time index
   # for dt = 1800 s (1/2 hour), there are 48 samples per day 
   # compute fsle for 5 days et =  start + 48 * 5
   et = st + (48 * 10) +1
   time = tr.time[st:et] - tr.time[st]    
   state = np.ones(np.size(tr.id))
   samples = tr.n_samples[:,st:et]
   criteria = 3

   latt = tr.lat[:,st:et]
   latt = np.squeeze(latt[:,timestep])
   latt = latt[np.where((state==1)&(~np.isnan(latt)))]

   lont = tr.lon[:,st:et]
   lont = np.squeeze(lont[:,timestep])
   lont = lont[np.where((state==1)&(~np.isnan(lont)))]

   uo = tr.u[:,st:et]
   uo = np.squeeze(uo[:,timestep])
   uo = uo[np.where((state==1)&(~np.isnan(uo)))]

   vo = tr.v[:,st:et]
   vo = np.squeeze(vo[:,timestep])
   vo = vo[np.where((state==1)&(~np.isnan(vo)))]

   lat0 = 28.8
   lon0 = -88.55
# Put observations in a 20 km by 20 km dimension
   xo = np.zeros(np.size(lont))
   yo = np.zeros(np.size(lont))
   for i in range(xo.size):
      xo[i] = GreatCircleDistance((latt[i],lont[i]),(latt[i],lon0)).km
      yo[i] = GreatCircleDistance((latt[i],lont[i]),(lat0,lont[i])).km
#      xo[i] = vincenty((latt[i],lont[i]),(latt[i],lon0)).km
#      yo[i] = vincenty((latt[i],lont[i]),(lat0,lont[i])).km
   obs = np.concatenate([uo,vo])
   obs = np.reshape(obs,[obs.size,1])

   dx = 0.5
   x = np.arange(0,20+dx,dx)
   y = np.arange(0,20+dx,dx)
   X,Y = np.meshgrid(x,y)
   Xs = np.reshape(X,[X.size])
   Ys = np.reshape(Y,[Y.size])   
#
   K = compute_K(xo,yo,sigma,1)
   Ki = np.linalg.inv(K)
   Ks = compute_Ks(xo,yo,Xs,Ys,sigma,1)
   f = getMean(Ks,Ki,obs)
   uf = np.reshape(f[:f.size/2],[y.size,-1])
   vf = np.reshape(f[f.size/2:],[y.size,-1])

   return x,y,uf,vf,xo,yo,uo,vo

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
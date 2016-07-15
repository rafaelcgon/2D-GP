import numpy as np
from datetime import datetime
from laser_class import drifter,interpolated_tracks
from laser_io_methods import save_object
import cPickle as pickle
import os
from dateutil.parser import parse
from datetime import datetime,timedelta
from subprocess import Popen, PIPE
from scipy import interpolate
import GPy
# open files
filename = 'ALL_2016_2_7.pkl'
kriging_filename = 'new_kriging2_' + filename
with open(kriging_filename,'rb') as input:
          kr = pickle.load(input)

N = np.size(kr.id)
u = kr.u
v = kr.v
drogStat = kr.drogueStat[:,:-1]
u[np.where(drogStat==0)]=np.nan
v[np.where(drogStat==0)]=np.nan
timem = (kr.time[1:]+kr.time[:-1])/2.

u2 = np.zeros((N,timem.size))
v2 = np.zeros((N,timem.size))

X = timem[::10]
Y1 = u[:,::10]
Y2 = v[:,::10]

kernel = GPy.kern.RBF(input_dim=1, variance=1159.68, lengthscale=4.5)
for n in range(N): 
    print n
    it = np.where(~np.isnan(Y1[n,:]))
    it2 = np.where((timem>=X[it[0][0]])&(timem<=X[it[0][-1]]))
    Y = np.squeeze(np.array([Y1[n,it],Y2[n,it]]).T)
    model = GPy.models.GPRegression(X[it],Y,kernel)
#    model.Gaussian_noise = 1.75598244486e-07 # got from previous experiments
    model.optimize_restarts()
    print model
    variables,vari = model.predict(timem[it2][:,None])
    u2[n,it2] = variables[:,0]
    v2[n,it2] = variables[:,1]   


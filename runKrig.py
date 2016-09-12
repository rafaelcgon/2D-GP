import numpy as np
import krig
import sys
import os
print sys.argv
print sys.argv[1]
ind = int(sys.argv[1])-1
print ind
T = np.array([1,2,2,3])
dt = np.array([2,2,4,4])
skp = 3
st = 0
et = T[ind]*24*60/15
outFile = 'rbfModel_T'+str(T[ind])+'_dt'+str(dt[ind])
outDir = '/scratch/projects/iskandarani/rafael/GP/skip_' + str(skp)
if not(os.path.isdir(outDir)):
   os.system('mkdir ' + outDir)
outFile = outDir+'/'+outFile

print outFile
krig.kriging(st,et,sample_step=-dt[ind],skip=skp,num_restarts=10,save=1,output=outFile)



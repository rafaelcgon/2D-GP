import numpy as np
import krig
import sys
print sys.argv
print sys.argv[1]
ind = int(sys.argv[1])-1
print ind


T = np.array([1,1,2,2,2,3,3,3])
dt = np.array([2,4,4,8,16,4,8,16])
skp = 5
st = 0
et = T[ind]*24*60/15
outFile = 'rbfModel_T'+str(T[ind])+'_dt'+str(dt[ind])
outDir = '/scratch/projects/iskandarani/rafael/GP/skip_' + str(skp)
outFile = outDir+'/'+outFile

krig.predict(outFile)

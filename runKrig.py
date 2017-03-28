import numpy as np
import krig
import sys
import os
print sys.argv
print sys.argv[1]
ind = int(sys.argv[1])-1
print ind
#T  =  np.array([1,2,3,4]) #2,3,3,4,4,5,20]) #1,1,2,2,3,3,4])
#dt  = np.array([3,4,5,5])#4,4,5,5,4,6,24]) #[1,2,2,4,2,4,4])
#skp = np.array([1,2,3,3,4])#3,3,3,3,4,3,3]) 
#nK = 2 # number of kernels (k1 + k2 +...+ Kn)

T  =  np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]) 
dt  = np.array([1,1,2,2,2,1,1,2,2,2,1,1,2,2,2]) 
skp = np.array([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3])  
nK =  np.array([1,2,1,2,3,1,2,1,2,3,1,2,1,2,3]) 

st = 0
et = st + T[ind]*24*60/15
# CHECK THE VALUE OF laser!!!!
#   1 for the data
#   0 for simulations
laser =0

outFile = 'rbfModel_T'+str(T[ind])+'_dt'+str(dt[ind])+'_nK'+str(nK[ind])
if laser==1:
   outDir = 'skip_' + str(skp[ind])
else:
   outDir = 'Simulations/skip_' + str(skp[ind])
if not(os.path.isdir(outDir)):
   os.system('mkdir ' + outDir)
outFile = outDir+'/'+outFile

print outFile
krig.kriging(st,et,sample_step=-dt[ind],skip=skp[ind],nKernels=nK[ind],output=outFile,laser=laser) #,pkg='scikit')



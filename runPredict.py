import numpy as np
import krig
import sys
print sys.argv
print sys.argv[1]
ind = int(sys.argv[1])-1
print ind


#T = np.array([1,1,2,2,2,3,3,3])
#dt = np.array([2,4,4,8,16,4,8,16])
#skp = 5
#st = 0
#et = T[ind]*24*60/15
#outFile = 'rbfModel_T'+str(T[ind])+'_dt'+str(dt[ind])
#outDir = '/scratch/projects/iskandarani/rafael/GP/skip_' + str(skp)
#outFile = outDir+'/'+outFile

#krig.predict(outFile)

T  = 1.0
dt = 2
skp= 2
nK = 2
#idt = np.array([-3,-2,-1,0,1,2,3,8,9,11])
idt0 = np.arange(8,24,4)
# idt = 12 
dtg = 0.5
tg = np.arange(12,36,dtg)
if ind>=idt0.size:
   ind = ind-idt0.size
   var='u'
else:
   var='v'
print 'Part ',ind+1,' of ',idt0.size,'.'

idt = idt0[ind]
simul = 0

ylim = [1,15]
xlim = [-5,15]
dx=0.1


outFile = 'rbfModel_T'+str(T) + '_dt'+str(dt)  + '_nK'+str(nK) 
outDir = 'outputs_proposal/skip2_day1/'
outFile = outDir+'/'+outFile
krig.scikit_prior(outFile,varname=var,dt=tg[idt],tlim=8,radar='',xlim=xlim,ylim=ylim,dx=dx)


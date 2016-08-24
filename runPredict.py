
import krig

T = np.array([1,2,2,3])
dt = np.array([2,2,4,4])
skp = 3
st = 0
et = T[ind]*24*60/15
outFile = 'rbfModel_T'+str(T[ind])+'_dt'+str(dt[ind])
outDir = '/scratch/projects/iskandarani/rafael/GP/skip_' + str(skp)
outFile = outDir+'/'+outFile

krig.predict(outFile)

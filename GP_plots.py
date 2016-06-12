import numpy as np
import scipy as sc
from scipy.interpolate import Rbf
from matplotlib import pylab as pl
import matplotlib.pyplot as plt               #
from matplotlib import rc,rcParams
from GP_scripts import *

rc('text',usetex=True)
####################################################################
def unit_vector(vector,ax=-1):
    """ Returns the unit vector of the vector.  """
    if ax == -1:
       return vector / np.linalg.norm(vector)
    else:
       nor = np.linalg.norm(vector,axis=1)
       nor = np.dot(np.reshape(nor,[nor.size,1]),np.array([[1,1]]))
       return vector / nor
##########################
def get_angle(v,ax=1):
    """ Returns the angle in degrees of a vector 'v'  ::

           - Vector v should have 2 components v[0]--> x comp; v[1]--> y comp.
           - For several vectors: 
                * the components should be arranged as vel = np.array([um,vm]).T
                  where um and vm are 1D arrays with the each vector component.
           - The northward unit vector is used as reference.
    """
    vu = unit_vector(v,ax)
    vn = np.array([[0,1]]).T 
    angle = np.degrees(np.arccos(np.clip(np.dot(vu, vn), -1.0, 1.0)))
    angle[np.where(v[:,0]<0)] = 360-angle[np.where(v[:,0]<0)] 

    return angle
####################################################################
def run_example(divFree = 1):
   x,y,phi,xm,ym,um,vm = generate_2D_gaussian(divFree)
   
   X1,X2=np.meshgrid(xm,ym)
   X1 = np.reshape(X1,[-1])
   X2 = np.reshape(X2,[-1])
   um = np.reshape(um,[-1])
   vm = np.reshape(vm,[-1])
# generate random samples
   samples = np.random.randint(0,X1.size,20)
   x1 = X1[samples]
   x2 = X2[samples]
   y1 = um[samples]
   y2 = vm[samples]
   y = np.concatenate([y1,y2])
   y = np.reshape(y,[y.size,1])

   x1s = xm #np.arange(x1.min(),x1.max(),0.1)
   x2s = ym #np.arange(x2.min(),x2.max(),0.1)
   X1s,X2s=np.meshgrid(x1s,x2s)   
   X1s = np.reshape(X1s,[X1s.size])
   X2s = np.reshape(X2s,[X2s.size])

   K1 = compute_K(x1,x2,0.2,divFree)
   Ki1 = np.linalg.inv(K1)
   KS1 = compute_Ks(x1,x2,X1s,X2s,0.2,divFree)
   f1 = getMean(KS1,Ki1,y)

   K2 = compute_K(x1,x2,0.3,divFree)
   Ki2 = np.linalg.inv(K2)
   KS2 = compute_Ks(x1,x2,X1s,X2s,0.3,divFree)
   f2 = getMean(KS2,Ki2,y)

   K3 = compute_K(x1,x2,0.4,divFree)
   Ki3 = np.linalg.inv(K3)
   KS3 = compute_Ks(x1,x2,X1s,X2s,0.4,divFree)
   f3 = getMean(KS3,Ki3,y)

   K4 = sqExp(x1,x2,x1,x2,0.3)
   Ki4 = np.linalg.inv(K4)
   KS4 = sqExp(X1s,X2s,x1,x2,0.3)
   u4 = np.dot(KS4,np.dot(Ki4,y1)) #getMean(KS4,Ki4,y1)
   u4 = np.reshape(u4,[x2s.size,-1])
   v4 = np.dot(KS4,np.dot(Ki4,y2)) #getMean(KS4,Ki4,y2)
   v4 = np.reshape(v4,[x2s.size,-1])

# errors by component
   rmsu1,rmsv1 = rmse(X1s,X2s,f1[:f1.size/2],f1[f1.size/2:],X1,X2,um,vm,knd = '')
   err_u1,ds = absoluteError(um,f1[:f1.size/2],X1s,X2s,x1,x2)
   err_v1,ds = absoluteError(vm,f1[f1.size/2:],X1s,X2s,x1,x2)

   rmsu2,rmsv2 = rmse(X1s,X2s,f2[:f2.size/2],f2[f2.size/2:],X1,X2,um,vm,knd = '')
   err_u2,ds = absoluteError(um,f2[:f2.size/2],X1s,X2s,x1,x2)
   err_v2,ds = absoluteError(vm,f2[f2.size/2:],X1s,X2s,x1,x2)

   rmsu3,rmsv3 = rmse(X1s,X2s,f3[:f3.size/2],f3[f3.size/2:],X1,X2,um,vm,knd = '')
   err_u3,ds = absoluteError(um,f3[:f3.size/2],X1s,X2s,x1,x2)
   err_v3,ds = absoluteError(vm,f3[f3.size/2:],X1s,X2s,x1,x2)

   rmsu4,rmsv4 = rmse(X1s,X2s,np.reshape(u4,[u4.size]),np.reshape(v4,[v4.size]),X1,X2,um,vm,knd = '')
   err_u4,ds = absoluteError(um,np.reshape(u4,[u4.size]),X1s,X2s,x1,x2)
   err_v4,ds = absoluteError(vm,np.reshape(v4,[v4.size]),X1s,X2s,x1,x2)

# errors by magnitude and angle

   absVel = np.sqrt(np.square(um)+np.square(vm))  
   angle =  get_angle(np.array([um,vm]).T)

   absVel1 = np.sqrt(np.square(f1[:f1.size/2])+np.square(f1[f1.size/2:]))   
   angle1 =  get_angle(np.array([f1[:f1.size/2],f1[f1.size/2:]]).T)
   err_vel1,ds = absoluteError(absVel,absVel1,X1s,X2s,x1,x2)
   err_ang1,ds = absoluteError(angle,angle1,X1s,X2s,x1,x2)

   absVel2 = np.sqrt(np.square(f2[:f2.size/2])+np.square(f2[f2.size/2:]))   
   angle2 =  get_angle(np.array([f2[:f2.size/2],f2[f2.size/2:]]).T)
   err_vel2,ds = absoluteError(absVel,absVel2,X1s,X2s,x1,x2)
   err_ang2,ds = absoluteError(angle,angle2,X1s,X2s,x1,x2)

   absVel3 = np.sqrt(np.square(f3[:f3.size/2])+np.square(f3[f3.size/2:]))   
   angle3 =  get_angle(np.array([f3[:f3.size/2],f3[f3.size/2:]]).T)
   err_vel3,ds = absoluteError(absVel,absVel3,X1s,X2s,x1,x2)
   err_ang3,ds = absoluteError(angle,angle3,X1s,X2s,x1,x2)

   absVel4 = np.sqrt(np.square(np.reshape(u4,[u4.size]))+np.square(np.reshape(v4,[v4.size])))   
   angle4 =  get_angle(np.array([np.reshape(u4,[u4.size]),np.reshape(v4,[v4.size])]).T)
   err_vel4,ds = absoluteError(absVel,absVel4,X1s,X2s,x1,x2)
   err_ang4,ds = absoluteError(angle,angle4,X1s,X2s,x1,x2)

   figW = 20.
   figH = 20.
   FS = 42
   FS2 = 20
   fig = pl.figure(figsize=(figW,figH))
   plot = fig.add_subplot(221)
   plot.plot(ds,err_u1,'ob',ds,err_u2,'or',ds,err_u3,'og',ds,err_u4,'ok')
   plt.legend(('div-free, $\sigma=0.1$','div-free, $\sigma=0.3$', 
     'div-free, $\sigma=0.4$','squared exp., $\sigma=0.3$'),
      fontsize=FS2,loc=1)
   plot.tick_params(axis='both',which='major',labelsize=FS2)
   plot.set_title('Zonal component',fontsize = FS)
   plot.set_ylabel('Absolute error',fontsize=FS)
   #plot.set_xlabel('Distance to nearest observation',fontsize=FS)

   plot = fig.add_subplot(222)
   plot.plot(ds,err_v1,'ob',ds,err_v2,'or',ds,err_v3,'og',ds,err_v4,'ok')
#   plt.legend(('div-free covariance, $\sigma=0.1$','div-free covariance, $\sigma=0.3$', 
#     'div-free covariance, $\sigma=0.4$','squared exponential covariance, $\sigma=0.3$'),
#      fontsize=FS2,loc=2)
   plot.tick_params(axis='both',which='major',labelsize=FS2)
   plot.set_title('Meridional component',fontsize = FS)
#   plot.set_ylabel('Absolute error',fontsize=FS)
#   plot.set_xlabel('Distance to nearest observation',fontsize=FS)

   plot = fig.add_subplot(223)
   plot.plot(ds,err_vel1,'ob',ds,err_vel2,'or',ds,err_vel3,'og',ds,err_vel4,'ok')
#   plt.legend(('div-free covariance, $\sigma=0.1$','div-free covariance, $\sigma=0.3$', 
#     'div-free covariance, $\sigma=0.4$','squared exponential covariance, $\sigma=0.3$'),
#      fontsize=FS2,loc=2)
   plot.tick_params(axis='both',which='major',labelsize=FS2)
   plot.set_title('Absolute velocity',fontsize = FS)
   plot.set_ylabel('Absolute error',fontsize=FS)
   plot.set_xlabel('Distance to nearest observation',fontsize=FS)

   plot = fig.add_subplot(224)
   plot.plot(ds,err_ang1,'ob',ds,err_ang2,'or',ds,err_ang3,'og',ds,err_ang4,'ok')
#   plt.legend(('div-free, $\sigma=0.1$','div-free, $\sigma=0.3$', 
#     'div-free, $\sigma=0.4$','squared exp., $\sigma=0.3$'),
#      fontsize=FS2,loc=4)
   plot.tick_params(axis='both',which='major',labelsize=FS2)
   plot.set_title('Angle',fontsize = FS)
#   plot.set_ylabel('Absolute error',fontsize=FS)
   plot.set_xlabel('Distance to nearest observation',fontsize=FS)

#   plot=fig.add_subplot(221)
#   plot.quiver(X1,X2,um,vm,scale = 1)
#   plot.quiver(x1,x2,y1,y2,scale=1,color='r')
#   plot.quiver(X1s,X2s,f1[0:f1.size/2],f1[f1.size/2:],scale = 1, color = 'b')
#   plot.set_title('$\sigma$ = 0.1',fontsize = FS)
#   plot.tick_params(axis='both',which='major',labelsize=FS2)

#   plot2=fig.add_subplot(222)
#   plot2.quiver(X1,X2,um,vm,scale = 1)
#   plot2.quiver(x1,x2,y1,y2,scale=1,color='r')
#   plot2.quiver(X1s,X2s,f2[0:f2.size/2],f2[f2.size/2:],scale = 1, color = 'b')
#   plot2.set_title('$\sigma$ = 0.3',fontsize = FS)
#   plot2.tick_params(axis='both',which='major',labelsize=FS2)

#   plot3=fig.add_subplot(223)
#   plot3.quiver(X1,X2,um,vm,scale = 1)
#   plot3.quiver(x1,x2,y1,y2,scale=1,color='r')
#   plot3.quiver(X1s,X2s,f3[0:f1.size/2],f3[f1.size/2:],scale = 1, color = 'b')
#   plot3.set_title('$\sigma$ = 0.4',fontsize = FS)
#   plot3.tick_params(axis='both',which='major',labelsize=FS2)

#   plot4=fig.add_subplot(224)
#   plot4.quiver(X1,X2,um,vm,scale = 1)
#   plot4.quiver(x1,x2,y1,y2,scale=1,color='r')
#   plot4.quiver(X1s,X2s,u4,v4,scale = 1, color = 'b')
#   plot4.set_title('Gaussian kernel, $\sigma$ = 0.3',fontsize = FS)
#   plot4.tick_params(axis='both',which='major',labelsize=FS2)

#   return X1,X2,um,vm,x1,x2,y1,y2,X1s,X2s,u4,v4
#####################################################################
#def plotAbsErr()

##############################3
def plotCovariance(sigma = 0.2,divFree=1,ds = 0.08):

   x = np.arange(-1,1+ds,ds)
   y = np.arange(-1,1+ds,ds)
   X,Y = np.meshgrid(x,y)
   x2 = np.reshape(X,[-1])
   y2 = np.reshape(Y,[-1])

   ori = x2.size/2 + x.size/2

   K = compute_K(x2,y2,sigma,divFree)
   Kxx = K[ori,:x2.size]
   Kxx = np.reshape(Kxx,[x.size,-1])
   Kxy = K[ori,x2.size:]
   Kxy = np.reshape(Kxy,[x.size,-1])
   Kyy = K[x2.size+ori,x2.size:]
   Kyy = np.reshape(Kyy,[x.size,-1])
 
   dxx = x-x2[ori]
   dyy = y-y2[ori]


   figW = 30.
   figH = 6.
   FS = 42
   FS2 = 20
   fig = pl.figure(figsize=(figW,figH))

   plot=fig.add_subplot(131)
   cs = plot.contourf(dxx,dxx,Kxx) 
   cbar = fig.colorbar(cs)
   plot.grid(); 
   plot.set_title('uu',fontsize = FS)
   plot.tick_params(axis='both',which='major',labelsize=FS2)
#   plot.axis('equal')
   plot.set_xlim([-1,1])
   plot.set_ylim([-1,1])

   plot=fig.add_subplot(132)
   cs = plot.contourf(dxx,dyy,Kxy) 
   cbar = fig.colorbar(cs)
   plot.grid(); 
   plot.set_title('uv',fontsize = FS)
   plot.tick_params(axis='both',which='major',labelsize=FS2)
#   plot.axis('equal')
   plot.set_xlim([-1,1])
   plot.set_ylim([-1,1])

   plot=fig.add_subplot(133)
   cs = plot.contourf(dyy,dyy,Kyy) 
   cbar = fig.colorbar(cs)
   plot.grid(); 
   plot.set_title('vv',fontsize = FS)
   plot.tick_params(axis='both',which='major',labelsize=FS2)
#   plot.axis('equal')
   plot.set_xlim([-1,1])
   plot.set_ylim([-1,1])
###############################################################################
def singleVector():
   dx = 0.05
   x = np.arange(-1,1+dx,dx)
   y = np.arange(-1,1+dx,dx)
   X,Y = np.meshgrid(x,y)
   xo = np.array([0.])
   yo = np.array([0.])
   uo = np.array([1.])
   vo = np.array([0.])
   ob = np.concatenate([uo,vo])
   ob = np.reshape(ob,[ob.size,1])
   K1 = gp.compute_K(xo,xo,0.1,1)
   Ks1 = gp.compute_Ks(xo,xo,X,Y,0.1,1)
   Ki1 = np.linalg.inv(K1)
   f1 = gp.getMean(Ks1,Ki1,ob)
   u1 = np.reshape(f1[0:f1.size/2],[x.size,-1])
   v1 = np.reshape(f1[f1.size/2:],[x.size,-1])










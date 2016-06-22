import numpy as np
import scipy as sc
from scipy.interpolate import Rbf
from matplotlib import pylab as pl
import matplotlib.pyplot as plt               #
from matplotlib import rc,rcParams
import matplotlib as mpl 
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
def plot4var(var1,var2,var3,var4,ds,title,text1,text2):
   maxvar = np.max([var1.max(),var2.max(),var3.max(),var4.max()])
   minvar = np.min([var1.min(),var2.min(),var3.min(),var4.min()])
   figW = 20.
   figH = 20.
   FS = 42
   FS2 = 35
   fig = pl.figure(figsize=(figW,figH))
   plot = fig.add_subplot(221)
   plot.plot(ds,var1,'ob')
   plot.set_ylim([minvar,maxvar])
   titleText = plot.text(0.8, 1.1, '', transform=plot.transAxes,fontsize=FS,fontweight='bold')
   titleText1 = plot.text(0.05, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText2 = plot.text(0.6, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText.set_text(title)
   titleText1.set_text(text1[0])
   titleText2.set_text(text2[0])
   plot.tick_params(axis='both',which='major',labelsize=FS2)
#   plot.set_title(title,fontsize = FS)
   plot.set_ylabel('Absolute error',fontsize=FS)
   #plot.set_xlabel('Distance to nearest observation',fontsize=FS)

   plot = fig.add_subplot(222)
   plot.plot(ds,var2,'or')
   plot.set_ylim([minvar,maxvar])
   titleText1 = plot.text(0.05, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText2 = plot.text(0.6, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText1.set_text(text1[1])
   titleText2.set_text(text2[1])
   plot.tick_params(axis='both',which='major',labelsize=FS2)
#   plot.set_title(title,fontsize = FS)

   plot = fig.add_subplot(223)
   plot.plot(ds,var3,'og')
   plot.set_ylim([minvar,maxvar])
   titleText1 = plot.text(0.05, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText2 = plot.text(0.6, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText1.set_text(text1[2])
   titleText2.set_text(text2[2])
   plot.tick_params(axis='both',which='major',labelsize=FS2)
#   plot.set_title(title,fontsize = FS)
   plot.set_ylabel('Absolute error',fontsize=FS)
   plot.set_xlabel('Distance to nearest obs.',fontsize=FS)

   plot = fig.add_subplot(224)
   plot.plot(ds,var4,'ok')
   plot.set_ylim([minvar,maxvar])
   titleText1 = plot.text(0.05, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText2 = plot.text(0.6, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText1.set_text(text1[3])
   titleText2.set_text(text2[3])
   plot.tick_params(axis='both',which='major',labelsize=FS2)
#   plot.set_title(title,fontsize = FS)
   plot.set_xlabel('Distance to nearest obs.',fontsize=FS)

   pl.savefig('plots/'+title+'-scatter.png', bbox_inches=0)
####################################################################
def contour4var(var1,var2,var3,var4,x,y,xo,yo,title,text1,text2):
   maxvar = np.max([var1.max(),var2.max(),var3.max(),var4.max()])
   minvar = np.min([var1.min(),var2.min(),var3.min(),var4.min()])
   dv = (maxvar-minvar)/32.
   V = np.arange(minvar,maxvar+dv,dv)

   var1 = np.reshape(var1,[x.size,-1])
   var2 = np.reshape(var2,[x.size,-1])
   var3 = np.reshape(var3,[x.size,-1])
   var4 = np.reshape(var4,[x.size,-1])

   figW = 20.
   figH = 14.
   bottom = [0.51,0.04]
   height = 0.4
   left = [0.04,0.5]
   width = 0.4

   FS = 42
   FS2 = 35
   fig = plt.figure(figsize=(figW,figH))
#   plot = fig.add_subplot(221)
   plot=fig.add_axes([left[0],bottom[0],width,height])
   cs=plot.contourf(x,y,var1,V)
   plot.plot(xo,yo,'ow',ms = 12)
   titleText = plot.text(0.8, 1.1, '', transform=plot.transAxes,fontsize=FS,fontweight='bold')
   titleText1 = plot.text(0.05, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText2 = plot.text(0.6, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText.set_text(title)
   titleText1.set_text(text1[0])
   titleText2.set_text(text2[0])
   plot.tick_params(axis='both',which='major',labelsize=FS2)
#   plot.set_title(title,fontsize = FS)

#   plot = fig.add_subplot(222)
   plot=fig.add_axes([left[1],bottom[0],width,height])
   cs=plot.contourf(x,y,var2,V)
   plot.plot(xo,yo,'ow',ms = 12)
   titleText1 = plot.text(0.05, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText2 = plot.text(0.6, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText1.set_text(text1[1])
   titleText2.set_text(text2[1])
   plot.tick_params(axis='both',which='major',labelsize=FS2)
#   plot.set_title(title,fontsize = FS)

#   plot = fig.add_subplot(223)
   plot=fig.add_axes([left[0],bottom[1],width,height])
   cs=plot.contourf(x,y,var3,V)
   plot.plot(xo,yo,'ow',ms = 12)
   titleText1 = plot.text(0.05, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText2 = plot.text(0.6, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText1.set_text(text1[2])
   titleText2.set_text(text2[2])
   plot.tick_params(axis='both',which='major',labelsize=FS2)
#   plot.set_title(title,fontsize = FS)

#   plot = fig.add_subplot(224)
   plot=fig.add_axes([left[1],bottom[1],width,height])
   cs = plot.contourf(x,y,var4,V)
   plot.plot(xo,yo,'ow',ms = 12)
   titleText1 = plot.text(0.05, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText2 = plot.text(0.6, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText1.set_text(text1[3])
   titleText2.set_text(text2[3])
   plot.tick_params(axis='both',which='major',labelsize=FS2)
#   plot.set_title(title,fontsize = FS)

## INDEPENDENT COLORBAR
   axbar = fig.add_axes([0.91, 0.1, 0.02, 0.8])
   cb = mpl.colorbar.ColorbarBase(axbar, orientation = 'vertical', boundaries = V)
   Vs = np.round(V[::4],decimals=3)
   cb.set_ticks(Vs)
   cb.ax.tick_params(labelsize = FS2)

   pl.savefig('plots/'+title+'-contour.png', bbox_inches=0)
#########################################################################
def quiver4var(u1,v1,u2,v2,u3,v3,u4,v4,x,y,xo,yo,title,text1,text2):

#   var1 = np.reshape(var1,[x.size,-1])
#   var2 = np.reshape(var2,[x.size,-1])
#   var3 = np.reshape(var3,[x.size,-1])
#   var4 = np.reshape(var4,[x.size,-1])

   FS = 42
   FS2 = 35
   scal = 1
   figW = 20.
   figH = 19.
   bottom = [0.51,0.04]
   height = 0.41
   left = [0.06,0.54]
   width = 0.41

   fig = pl.figure(figsize=(figW,figH))
#   plot = fig.add_subplot(221)
   plot=fig.add_axes([left[0],bottom[0],width,height])
   cs=plot.quiver(x,y,u1,v1,scale=scal)
   plot.plot(xo,yo,'or',ms = 15)
   titleText = plot.text(0.85, 1.1, '', transform=plot.transAxes,fontsize=FS,fontweight='bold')
   titleText1 = plot.text(0.05, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText2 = plot.text(0.6, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText.set_text(title)
   titleText1.set_text(text1[0])
   titleText2.set_text(text2[0])
   plot.tick_params(axis='both',which='major',labelsize=FS2)
#   plot.set_title(title,fontsize = FS)

#   plot = fig.add_subplot(222)
   plot=fig.add_axes([left[1],bottom[0],width,height])
   cs=plot.quiver(x,y,u2,v2,scale=scal)
   plot.plot(xo,yo,'or',ms = 15)
   titleText1 = plot.text(0.05, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText2 = plot.text(0.6, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText1.set_text(text1[1])
   titleText2.set_text(text2[1])
   plot.tick_params(axis='both',which='major',labelsize=FS2)
#   plot.set_title(title,fontsize = FS)

#   plot = fig.add_subplot(223)
   plot=fig.add_axes([left[0],bottom[1],width,height])
   cs=plot.quiver(x,y,u3,v3,scale=scal)
   plot.plot(xo,yo,'or',ms = 15)
   titleText1 = plot.text(0.05, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText2 = plot.text(0.6, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText1.set_text(text1[2])
   titleText2.set_text(text2[2])
   plot.tick_params(axis='both',which='major',labelsize=FS2)
#   plot.set_title(title,fontsize = FS)

#   plot = fig.add_subplot(224)
   plot=fig.add_axes([left[1],bottom[1],width,height])
   cs=plot.quiver(x,y,u4,v4,scale=scal)
   plot.plot(xo,yo,'or',ms = 15)
   titleText1 = plot.text(0.05, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText2 = plot.text(0.6, 1.01, '', transform=plot.transAxes,fontsize=FS2,fontweight='bold',color='k')
   titleText1.set_text(text1[3])
   titleText2.set_text(text2[3])
   plot.tick_params(axis='both',which='major',labelsize=FS2)
#   plot.set_title(title,fontsize = FS)

   pl.savefig('plots/'+title+'.png', bbox_inches=0)
#####################################################################
def run_example(nsamples = 5,divFree = 1):
   x,y,phi,xm,ym,um,vm = generate_2D_gaussian(divFree)
   
   X1,X2=np.meshgrid(xm,ym)
   X1 = np.reshape(X1,[-1])
   X2 = np.reshape(X2,[-1])
   um = np.reshape(um,[-1])
   vm = np.reshape(vm,[-1])
# generate random samples
   
   samples = np.random.randint(0,X1.size,nsamples)
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
   
   sigma = np.array([0.1,0.3,0.5])
   K1 = compute_K(x1,x2,sigma[0],divFree)
   Ki1 = np.linalg.inv(K1)
   KS1 = compute_Ks(x1,x2,X1s,X2s,sigma[0],divFree)
   f1 = getMean(KS1,Ki1,y)
   u1 = np.reshape(f1[:f1.size/2],[x2s.size,-1])
   v1 = np.reshape(f1[f1.size/2:],[x2s.size,-1])

   K2 = compute_K(x1,x2,sigma[1],divFree)
   Ki2 = np.linalg.inv(K2)
   KS2 = compute_Ks(x1,x2,X1s,X2s,sigma[1],divFree)
   f2 = getMean(KS2,Ki2,y)
   u2 = np.reshape(f2[:f2.size/2],[x2s.size,-1])
   v2 = np.reshape(f2[f2.size/2:],[x2s.size,-1])

   K3 = compute_K(x1,x2,sigma[2],divFree)
   Ki3 = np.linalg.inv(K3)
   KS3 = compute_Ks(x1,x2,X1s,X2s,sigma[2],divFree)
   f3 = getMean(KS3,Ki3,y)
   u3 = np.reshape(f3[:f3.size/2],[x2s.size,-1])
   v3 = np.reshape(f3[f3.size/2:],[x2s.size,-1])

   K4 = sqExp(x1,x2,x1,x2,sigma[2])
   Ki4 = np.linalg.inv(K4)
   KS4 = sqExp(X1s,X2s,x1,x2,sigma[2])
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
   err_ang1[np.where(err_ang1>180)]=360-err_ang1[np.where(err_ang1>180)]

   absVel2 = np.sqrt(np.square(f2[:f2.size/2])+np.square(f2[f2.size/2:]))   
   angle2 =  get_angle(np.array([f2[:f2.size/2],f2[f2.size/2:]]).T)
   err_vel2,ds = absoluteError(absVel,absVel2,X1s,X2s,x1,x2)
   err_ang2,ds = absoluteError(angle,angle2,X1s,X2s,x1,x2)
   err_ang2[np.where(err_ang2>180)]=360-err_ang2[np.where(err_ang2>180)]

   absVel3 = np.sqrt(np.square(f3[:f3.size/2])+np.square(f3[f3.size/2:]))   
   angle3 =  get_angle(np.array([f3[:f3.size/2],f3[f3.size/2:]]).T)
   err_vel3,ds = absoluteError(absVel,absVel3,X1s,X2s,x1,x2)
   err_ang3,ds = absoluteError(angle,angle3,X1s,X2s,x1,x2)
   err_ang3[np.where(err_ang3>180)]=360-err_ang3[np.where(err_ang3>180)]

   absVel4 = np.sqrt(np.square(np.reshape(u4,[u4.size]))+np.square(np.reshape(v4,[v4.size])))   
   angle4 =  get_angle(np.array([np.reshape(u4,[u4.size]),np.reshape(v4,[v4.size])]).T)
   err_vel4,ds = absoluteError(absVel,absVel4,X1s,X2s,x1,x2)
   err_ang4,ds = absoluteError(angle,angle4,X1s,X2s,x1,x2)
   err_ang4[np.where(err_ang4>180)]=360-err_ang4[np.where(err_ang4>180)]

   text1 = np.array(['Div-free Kernel','Div-free kernel','Div-free kernel','Isotropic kernel'])
   text2 = np.array(['$\sigma='+str(sigma[0])+'$',
                     '$\sigma='+str(sigma[1])+'$',
                     '$\sigma='+str(sigma[2])+'$',
                     '$\sigma='+str(sigma[2])+'$'])
   
#   plot4var(err_u1,err_u2,err_u3,err_u4,ds,'Absolute-Error-U',text1,text2)
#   plot4var(err_v1,err_v2,err_v3,err_v4,ds,'Absolute-Error-V',text1,text2)
#   plot4var(err_vel1,err_vel2,err_vel3,err_vel4,ds,'Absolute-Error-Speed',text1,text2)
#   plot4var(err_ang1,err_ang2,err_ang3,err_ang4,ds,'Absolute-Error-Angle',text1,text2)

#   contour4var(err_u1,err_u2,err_u3,err_u4,x1s,x2s,x1,x2,'Absolute-Error-U',text1,text2)
#   contour4var(err_v1,err_v2,err_v3,err_v4,x1s,x2s,x1,x2,'Absolute-Error-V',text1,text2)
#   contour4var(err_vel1,err_vel2,err_vel3,err_vel4,x1s,x2s,x1,x2,'Absolute-Error-Speed',text1,text2)
#   contour4var(err_ang1,err_ang2,err_ang3,err_ang4,x1s,x2s,x1,x2,'Absolute-Error-Angle',text1,text2)
   
   um = np.reshape(um,[xm.size,-1])
   vm = np.reshape(vm,[xm.size,-1])
#   return x1s,x2s,u1,v1,u2,v2,u3,v3,u4,v4,um,vm,x1,x2 
#   quiver4var(u1-um,v1-vm,u2-um,v2-vm,u3-um,v3-vm,u4-um,v4-vm,X1s,X2s,x1,x2,'Velocity-Error',text1,text2)
   quiver4var(u1,v1,u2,v2,u3,v3,u4,v4,X1s,X2s,x1,x2,'Reconstructed-Velocity',text1,text2)

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
   x2 = np.reshape(X,[-1])
   y2 = np.reshape(Y,[-1])
   xo = np.array([0.])
   yo = np.array([0.])
   uo = np.array([1.])
   vo = np.array([0.])
   ob = np.concatenate([uo,vo])
   ob = np.reshape(ob,[ob.size,1])
   
   K1 = sqExp(xo,xo,xo,xo,0.1)
   Ki1 = np.linalg.inv(K1)
   Ks1 = sqExp(x2,y2,xo,yo,0.1)
   u1 = np.dot(Ks1,np.dot(Ki1,uo)) 
   u1 = np.reshape(u1,[x.size,-1])
   v1 = np.dot(Ks1,np.dot(Ki1,vo)) 
   v1 = np.reshape(v1,[x.size,-1])

   K2 = compute_K(xo,xo,0.1,1)
   Ks2 = compute_Ks(xo,xo,x2,y2,0.1,1)
   Ki2 = np.linalg.inv(K2)
   f2 = getMean(Ks2,Ki2,ob)
   u2 = np.reshape(f2[0:f2.size/2],[x.size,-1])
   v2 = np.reshape(f2[f2.size/2:],[x.size,-1])

   K3 = compute_K(xo,xo,0.1,2)
   Ks3 = compute_Ks(xo,xo,x2,y2,0.1,2)
   Ki3 = np.linalg.inv(K3)
   f3 = getMean(Ks3,Ki3,ob)
   u3 = np.reshape(f3[0:f3.size/2],[x.size,-1])
   v3 = np.reshape(f3[f3.size/2:],[x.size,-1])

   figH = 20
   figW = 9.5
   bottom = [0.69,0.36,0.03]
   height = 0.27
   left = 0.085
   width = 0.88

   FS = 38
   FS2 = 32
   scal = 15.

   fig = pl.figure(figsize=(figW,figH))

   plot=fig.add_axes([left,bottom[0],width,height])
   cs=plot.quiver(x,y,u1,v1,scale=scal)
   pl.quiver(xo,yo,uo,vo,scale=scal,color='r')
   plot.tick_params(axis='both',which='major',labelsize=FS2)
   plot.set_title('Isotropic kernel',fontsize = FS)
   plot.set_xlim([-0.3,0.3])
   plot.set_ylim([-0.3,0.3])
   plot.set_xticks([-0.3,-0.2,-0.1,0.,0.1,0.2,0.3])
   plot.set_yticks([-0.2,-0.1,0.0,0.1,0.2])

   plot=fig.add_axes([left,bottom[1],width,height])
   cs=plot.quiver(x,y,u2,v2,scale=scal)
   pl.quiver(xo,yo,uo,vo,scale=scal,color='r')
   plot.tick_params(axis='both',which='major',labelsize=FS2)
   plot.set_title('Divergence-free kernel',fontsize = FS)
   plot.set_xlim([-0.3,0.3])
   plot.set_ylim([-0.3,0.3])
   plot.set_xticks([-0.3,-0.2,-0.1,0.,0.1,0.2,0.3])
   plot.set_yticks([-0.2,-0.1,0.0,0.1,0.2])

   plot=fig.add_axes([left,bottom[2],width,height])
   cs=plot.quiver(x,y,u3,v3,scale=scal)
   pl.quiver(xo,yo,uo,vo,scale=scal,color='r')
   plot.tick_params(axis='both',which='major',labelsize=FS2)
   plot.set_title('Curl-free kernel',fontsize = FS)
   plot.set_xlim([-0.3,0.3])
   plot.set_ylim([-0.3,0.3])
   plot.set_xticks([-0.3,-0.2,-0.1,0.,0.1,0.2,0.3])
   plot.set_yticks([-0.2,-0.1,0.0,0.1,0.2])

   pl.savefig('plots/one_obs.png', bbox_inches=0)
###############################################################################
def twoVectors(sigma = 0.2):
#
   figH = 18
   figW = 13
   bottom = [0.68,0.35,0.02]
   height = 0.3
   left = [0.07,0.56]
   width = 0.4

   FS = 42
   FS2 = 32

   scal = 25.
#
   dx = 0.05
   x = np.arange(-1,1+dx,dx)
   y = np.arange(-1,1+dx,dx)
   X,Y = np.meshgrid(x,y)
   Xs = np.reshape(X,[X.size])
   Ys = np.reshape(Y,[Y.size])   
# observation 1, at origin
   xo1 = 0.
   yo1 = 0.
   uo1 = 1.
   vo1 = 0.
# observation 2, varying direction, intensity and location
#   xo2 = np.array([0.07,0.14,0.28])
#   yo2 = np.array([0.07,0.14,0.28])
   xo2 = np.array([0.07,0.28])
   yo2 = np.array([0.07,0.28])
   angle = np.array([45,90,180])
   absVel = np.array([0.5,1,2])
#
   for i in range(angle.size):
      fig = pl.figure(figsize=(figW,figH))
      jk=1
      filename = 'plots/2_vectors_angle_'+str(angle[i])+'.png'
      for j in range(absVel.size):
         uo2 = absVel[j]*np.cos(np.deg2rad(angle[i]))   
         vo2 = absVel[j]*np.sin(np.deg2rad(angle[i]))   
         uo = np.array([uo1,uo2])
         vo = np.array([vo1,vo2])
         ob = np.concatenate([uo,vo])
         ob = np.reshape(ob,[ob.size,1])
         for k in range(xo2.size):
            xo = np.array([xo1,xo2[k]])
            yo = np.array([yo1,yo2[k]])
            # compute kernel    
            K = compute_K(xo,yo,sigma)
            Ki = np.linalg.inv(K)
            Ks = compute_Ks(xo,yo,Xs,Ys,sigma)
            f = getMean(Ks,Ki,ob)
            uf = np.reshape(f[:f.size/2],[y.size,-1])
            vf = np.reshape(f[f.size/2:],[y.size,-1])
            # plot
            print jk
            print i,j,k
#            plot = fig.add_subplot(3,3,jk)
            plot=fig.add_axes([left[k],bottom[j],width,height])
            plot.quiver(X,Y,uf,vf,scale=scal)
            pl.quiver(xo,yo,uo,vo,scale=scal,color='r')
            plot.set_xlim([-0.6,0.8])
            plot.set_ylim([-0.6,0.8])
            plot.set_xticks([-0.6,-0.3,0.,0.3,0.6])
            plot.set_yticks([-0.5,-0.2,0.1,0.4,0.7])

#            titleText = plot.text(0.7, 1.05, '',transform=plot.transAxes,
#                               fontsize=FS,fontweight='bold')
            titleText1 = plot.text(0.05, 0.9, '',transform=plot.transAxes,
                               fontsize=FS2,fontweight='bold',color='b')
            titleText2 = plot.text(0.03, 0.03, '',transform=plot.transAxes,
                               fontsize=FS2,fontweight='bold',color='b')
            titleText3 = plot.text(0.7, 0.03, '',transform=plot.transAxes,
                               fontsize=FS2,fontweight='bold',color='b')
            titleText3.set_text('$\Theta='+str(angle[i])+'$')
            titleText2.set_text('$|v_2|='+str(absVel[j])+'$')
            dx2= np.round(np.sqrt(np.square(xo1-xo2[k])+np.square(yo1-yo2[k])),decimals=2)
            titleText1.set_text('$|x_2-x_1|='+str(dx2)+'$')
            plot.tick_params(axis='both',which='major',labelsize=FS2)
            jk+=1
      pl.savefig(filename, bbox_inches=0)

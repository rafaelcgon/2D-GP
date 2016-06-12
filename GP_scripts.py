import numpy as np
import scipy as sc
from scipy.interpolate import Rbf

#################################################################################
def getMean(KS,Ki,y):
   f = np.dot(KS,np.dot(Ki,y))
   return np.reshape(f,[-1])
#################################################################################
def getMargLikelihood(x1,x2,x1s,x2s,sigma=0.2,divFree=1):
   K = compute_K(x1,x2,sigma,divFree)
   Ki = np.linalg.inv(K)
   Ks = compute_Ks(x1,x2,x1s,x2s,sigma,divFree)
   Kss = compute_K(x1s,x2s,sigma,divFree)
   ML = Kss - np.dot(Ks,np.dot(Ki,Ks.T))
   return ML,Ki,Ks

#####################################################
def nonDivK(xa,xb,sigma,divFree=1):
   p = np.size(xa) # number of dimensions
   B = np.zeros((p,1)) # vector used to compute the A matrix
   C = np.square(np.linalg.norm(xa-xb)/sigma) # Exponential  
   for i in range(p):
      B[i] = (xa[i] - xb[i])/sigma
   if divFree==1: # divergence free
      A = B*B.T + (((p-1) - C)*np.identity(p))
   elif divFree==2: # curl free
      A = np.identity(p) - B*B.T
   else:
      A = np.square(sigma) # for a simple square exponential
   return np.square(1./sigma)*np.exp(-C/2.)*A # D by D K matrix 
#####################################################
def compute_K(x1,x2,sigma,divFree=1): # compute K matrix at sample inputs x1,x2
   N = np.size(x1) # number of samples
   D = 2  # dimensions
   K = np.zeros((N*D,N*D))   
   ii=0 
   jj=0
   for i in range(N):
 #     print i,'/',N
      for j in range(i,N):
         Kij = nonDivK(np.array([x1[i],x2[i]]),np.array([x1[j],x2[j]]),sigma,divFree) 
         K[ii:ii+D,jj:jj+D] = Kij
         jj+=2
      ii+=2
      jj=ii
      K[ii:ii+D,0:jj] = K[0:jj,ii:ii+D].T # cp upper off-diag. elem. into lower off-diag. 
# map K from N x N blocks to D x D 
   K2 = np.zeros((N*D,N*D))
   K2[0:N,0:N]     = K[0:N*D:D,0:N*D:D]
   K2[N:2*N,0:N]   = K[1:N*D:D,0:N*D:D]
   K2[0:N,N:2*N]   = K[0:N*D:D,1:N*D:D]
   K2[N:2*N,N:2*N] = K[1:N*D:D,1:N*D:D]
   return K2
##########################################################################
def compute_Ks(x1,x2,x1s,x2s,sigma,divFree=1):
   N = np.size(x1)
   M = np.size(x1s)
   # turn x1s and x2s into 1D arrays
   x1s=np.reshape(x1s,[-1])
   x2s=np.reshape(x2s,[-1])

   D = 2
   Ks = np.zeros((M*D,N*D)) # it's one Ks matrix for each grid cell  
   ii=0
   jj=0
   for i in range(M):
      for j in range(N):
         Kij = nonDivK(np.array([x1s[i],x2s[i]]),np.array([x1[j],x2[j]]),sigma,divFree) 
         #print Kij,ns,ii 
         Ks[ii:ii+D,jj:jj+D] = Kij
         jj+=2
      ii+=2
      jj=0

   K2 = np.zeros((M*D,N*D))
   K2[0:M,0:N]     = Ks[0:M*D:D,0:N*D:D]
   K2[M:D*M,0:N]   = Ks[1:M*D:D,0:N*D:D]
   K2[0:M,N:D*N]   = Ks[0:M*D:D,1:N*D:D]
   K2[M:D*M,N:D*N] = Ks[1:M*D:D,1:N*D:D]
  
   return K2 
################################################################################
def sqExp(x1,y1,x2,y2,sigma):
    I = x1.size
    J = x2.size
    K = np.zeros((I,J))
    for i in range(I):
        for j in range(J):
            xa = np.array([x1[i],y1[i]])
            xb = np.array([x2[j],y2[j]])
            K[i,j] = nonDivK(xa,xb,sigma,divFree=0)
    return K
#################################################################################
def rmse(x1s,x2s,f1,f2,x1,x2,y1,y2,knd = ''):
   # compute the rmse in the data points
   # first: interpolate f into sample points
   if (knd=='cubic')|(knd=='linear'):
      f1g = sc.interpolate.interp2d(x1s, x2s, f1, kind = knd)
      f2g = sc.interpolate.interp2d(x1s, x2s, f2, kind = knd)
      y1s = np.zeros(x1.size)
      y2s = np.zeros(x1.size)
#      for i in range(x1.size):
#         y1s[i]=f1g(x1[i],x2[i])
#         y2s[i]=f2g(x1[i],x2[i])
      y1s=f1g(x1,x2)
      y2s=f2g(x1,x2)
   elif knd == 'rbf':
      f1g = Rbf(x1s, x2s, f1)
      f2g = Rbf(x1s, x2s, f2)
      y1s = f1g(x1,x2)
      y2s = f2g(x1,x2)
   else:
      y1s = f1
      y2s = f2

   # compute error for each component
   error1 = y1s-y1
   error1 = np.reshape(error1,[error1.size])
   rmse1 = np.sqrt(np.mean(np.square(error1)))
   error2 = y2s-y2
   error2 = np.reshape(error2,[error2.size])
   rmse2 = np.sqrt(np.mean(np.square(error2)))
   print rmse1/(y1.max()-y1.min()),rmse2/(y2.max()-y2.min())
   return rmse1,rmse2
################################################################################
def absoluteError(y,f,x1f,x2f,x1,x2):
   # y is the real field (scalar), 
   # f is the reconstructed field in the same coord as y
   # x1f and x2f are vectors with the position of each f
   # x1 and x2 are the positions of the observations used to compute 
   # the distance of each grid point to the nearest data point
   f = np.reshape(f,[-1])
   y = np.reshape(y,[-1])
   abs_error = np.abs(y-f)
#   ds = np.zeros(x1f.size)
#   for i in range(x1f.size):
#      dsi = np.sqrt(np.square(x1 - x1f[i]) + np.square(x2 - x2f[i]))
#      ds[i]= dsi.min()
 
   X1f,X1 = np.meshgrid(x1f,x1)
   X2f,X2 = np.meshgrid(x2f,x2)
   dsi = np.sqrt(np.square(X1 - X1f) + np.square(X2 - X2f))
   ds = np.min(dsi,0)
   return abs_error,ds
#################################################################################
def generate_2D_gaussian(divFree=1):
    dx = 0.05
    x = np.arange(-1,1+dx,dx)
    dy = 0.05
    y = np.arange(-1,1+dy,dy)
    A = 1.5
    l = 3.
    X,Y = np.meshgrid(x,y)
    phi = A*np.exp(-(X**2)/l - (Y**2)/l)
    dphi_dy = np.diff(phi,axis=0)
    dphi_dx = np.diff(phi,axis=1)

    if (divFree==1):
       um = (dphi_dy[:,:-1]+dphi_dy[:,1:])/2.
       vm = -(dphi_dx[:-1,:]+dphi_dx[1:,:])/2.
    else:
       um = -(dphi_dx[:-1,:]+dphi_dx[1:,:])/2.
       vm = -(dphi_dy[:,:-1]+dphi_dy[:,1:])/2.
    xm = (x[:-1]+x[1:])/2.
    ym = (x[:-1]+x[1:])/2.
    
    return x,y,phi,xm,ym,um,vm

#####################################################################################
def vel_grad(x,y,u,v):
    dx = np.diff(x)
    dy = np.diff(y)
    du = np.diff(u,axis=1)
    dv = np.diff(v,axis=0)
    dum = (du[1:,:]+du[:-1,:])/2.
    dvm = (dv[:,1:]+dv[:,:-1])/2.
    return dum/dx + dvm/dy


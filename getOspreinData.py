from netCDF4 import Dataset
import numpy as np
from scipy import interpolate
import matplotlib as mpl
from matplotlib import pylab as pl
from matplotlib import animation
import matplotlib.pyplot as plt               #
from mpl_toolkits.basemap import Basemap
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Getting matrix index ######################################################
def getIndices(Z,Y,X,latlim,lonlim,zlim):
   IZ = range(np.min(np.where(Z>=zlim[0])),np.max(np.where(Z<=zlim[1]))+ 1); NZ = np.size(IZ);
#   print np.size(np.where(Y>=latlim[0])),np.size(np.where(Y<=latlim[1]))
#   print latlim,Y.max()
   IY = range(np.min(np.where(Y>=latlim[0])),np.max(np.where(Y<=latlim[1])) + 1); NY = np.size(IY);
   IX = range(np.min(np.where(X>=lonlim[0])),np.max(np.where(X<=lonlim[1])) + 1); NX = np.size(IX);
   IT = range(0, 10); NT = 10;
   return IT,IZ,IY,IX,NT,NZ,NY,NX
####################################################################################################
def getVel(inFileName,latlim,lonlim,zlim,NCOM=1):
   print "Extracting data from "+inFileName
   inFile = Dataset(inFileName, 'r')
   X = inFile.variables['lon'][:]
   Y = inFile.variables['lat'][:]
   Z = inFile.variables['depth'][:]
   T = inFile.variables['time'][:]
   IT,IZ,IY,IX,NT,NZ,NY,NX = getIndices(Z,Y,X,latlim,lonlim,zlim)
   if NCOM ==1:
      vvel = inFile.variables['vvel'][IT,IZ,IY,IX]
      uvel = inFile.variables['uvel'][IT,IZ,IY,IX]
   else:
      vvel = inFile.variables['vvel0'][IT,IZ,IY,IX]
      uvel = inFile.variables['uvel0'][IT,IZ,IY,IX]
   return Z[IZ],Y[IY],X[IX],T[IT],uvel,vvel
####################################################################################################
def getSal(inFileName,latlim,lonlim,zlim):
   print "Extracting data from "+inFileName
   inFile = Dataset(inFileName, 'r')
   X = inFile.variables['lon'][:]
   Y = inFile.variables['lat'][:]
   Z = inFile.variables['depth'][:]
   T = inFile.variables['time'][:]
   IT,IZ,IY,IX,NT,NZ,NY,NX = getIndices(Z,Y,X,latlim,lonlim,zlim)
   sal = inFile.variables['saln'][IT,IZ,IY,IX]
   return Z[IZ],Y[IY],X[IX],T[IT],sal

# NetCDF ###########################################################################

def getOsprein(period = [2,29],latlim=[23.0,31.0],lonlim=[-91.,-82.0],zlim=[0,0],NCOM=0):
    # period = days of the year (2 - 29 is 01/02/2013 - 01/29/2013)
   #print 'Extracting data:'
   if NCOM ==1:
      inFileName1 = '/home/rgoncalves/NCOM_1km/osprein/osprein_2013_'
   else:
      inFileName1 = '/home/rgoncalves/LagOil/OSPRE/dcom/data/eofdata/eof_2010_'

   for n in range(period[0],period[1]+1):
       if NCOM ==1:
          inFileName =inFileName1  + str(n)+ '.nc'
       else: 
          inFileName =inFileName1  + str(n+108)+ '.nc'

   #    print inFileName
       z,y,x,time1,uvel1,vvel1 = getVel(inFileName,latlim,lonlim,zlim,NCOM)
       uvel1[np.where(uvel1>10000)]=0.
       vvel1[np.where(vvel1>10000)]=0.

       NT = np.size(time1)
       if (n==period[0]):
          if (n<period[1]): #Save the last data of the previous period
             #print 'first data'
             uvel = uvel1[0:NT-1,:,:,:]
             vvel= vvel1[0:NT-1,:,:,:] 
             time= time1[0:NT-1]
          else:
             uvel = uvel1[0:NT,:,:,:]
             vvel= vvel1[0:NT,:,:,:] 
             time= time1[0:NT]
       elif(n < period[1]):
          uvel = np.concatenate((uvel,uvel1[1:NT-1,:,:,:]),axis=0)
          vvel = np.concatenate((vvel,vvel1[1:NT-1,:,:,:]),axis=0)
          time = np.concatenate((time,time1[1:NT-1]),axis=0)
       elif (n>0): #save the first data of the next period 
          uvel = np.concatenate((uvel,uvel1[1:NT,:,:,:]),axis=0)
          vvel = np.concatenate((vvel,vvel1[1:NT,:,:,:]),axis=0)
          time = np.concatenate((time,time1[1:NT]),axis=0)
       else:
          print 'should not get this message '
          uvel = uvel1[0:NT,:,:,:]
          vvel=vvel1[0:NT,:,:,:] 
          time=time1[0:NT]
   return uvel,vvel,time,x,y

###########################################################################################
def getSalinity(period = [2,29],latlim=[23.0,31.0],lonlim=[-91.,-82.0],zlim=[0,0]):
    # period = days of the year (2 - 29 is 01/02/2013 - 01/29/2013)
   #print 'Extracting data:'
   inFileName1 = '/home/rgoncalves/NCOM_1km/osprein_w_saln/osprein_2013_'
   for n in range(period[0],period[1]+1):
       inFileName =inFileName1  + str(n)+ '.nc'
   #    print inFileName
       z,y,x,time1,sal1 = getSal(inFileName,latlim,lonlim,zlim)
       sal1[np.where(sal1>10000)]=0.
       
       NT = np.size(time1)
       if (n==period[0]):
          if (n<period[1]): #Save the last data of the previous period
             #print 'first data'
             sal = sal1[0:NT-1,:,:,:]
             time= time1[0:NT-1]
          else:
             sal = sal1[0:NT,:,:,:]
             time= time1[0:NT]
       elif(n < period[1]):
          sal = np.concatenate((sal,sal1[1:NT-1,:,:,:]),axis=0)
          time = np.concatenate((time,time1[1:NT-1]),axis=0)
       elif (n>0): #save the first data of the next period 
          sal = np.concatenate((sal,sal1[1:NT,:,:,:]),axis=0)
          time = np.concatenate((time,time1[1:NT]),axis=0)
       else:
          print 'should not get this message '
          sal = sal1[0:NT,:,:,:]
          time=time1[0:NT]
   x,y = np.meshgrid(x,y)
   return np.squeeze(sal),x,y,time

###########################################################################################

def getSalGrad(period = [2,29],latlim=[23.0,31.0],lonlim=[-91.,-82.0],zlevel=0):
    sal,lons,lats,times = getSalinity(period,latlim,lonlim,[zlevel,zlevel])
    lat0,time0,lon0 = np.meshgrid(lats[:,0],times,lons[0,:])
    print 'Computing salinity gradient'
    x0 = (lon0-lons[0,0])*(40000./360.)*np.cos(np.deg2rad(lat0))
    y0 = (lat0-lats[0,0])*(40000./360.)

    dsdx=np.diff(np.squeeze(sal),axis=2)/np.diff(x0,axis=2)
    dsdx = (dsdx[:,0:-1,:]+dsdx[:,1:,:])/2.
    dsdy=np.diff(np.squeeze(sal),axis=1)/np.diff(y0,axis=1)
    dsdy = (dsdy[:,:,0:-1]+dsdy[:,:,1:])/2.
    gradSal = dsdx + dsdy
    lon0 = np.squeeze((lon0[0,:,0:-1]+lon0[0,:,1:])/2.)   
    lon0 = np.squeeze((lon0[0:-1,:]+lon0[1:,:])/2.)   
    lat0 = np.squeeze((lat0[0,:,0:-1]+lat0[0,:,1:])/2.)   
    lat0 = np.squeeze((lat0[0:-1,:]+lat0[1:,:])/2.)   
    time0= np.squeeze(time0[:,0,0])

    return gradSal,lon0,lat0,time0

###########################################################################################

def getDiv(period = [2,29],latlim=[23.0,31.0],lonlim=[-91.,-82.0],zlevel=0):

    uvel,vvel,time_v,lon_v,lat_v=getOsprein(period,latlim,lonlim,[zlevel,zlevel])
    lat0,time0,lon0 = np.meshgrid(lat_v,time_v,lon_v)
    print 'Computing velocity divergence'
    x0 = (lon0-lon_v[0])*(40000./360.)*np.cos(np.deg2rad(lat0))
    y0 = (lat0-lat_v[0])*(40000./360.)

    dudx=np.diff(np.squeeze(uvel),axis=2)/np.diff(x0,axis=2)
    dudx = (dudx[:,0:-1,:]+dudx[:,1:,:])/2.
    dvdy=np.diff(np.squeeze(vvel),axis=1)/np.diff(y0,axis=1)
    dvdy = (dvdy[:,:,0:-1]+dvdy[:,:,1:])/2.
    velDiv0 = dudx + dvdy
    lon0 = np.squeeze((lon0[0,:,0:-1]+lon0[0,:,1:])/2.)   
    lon0 = np.squeeze((lon0[0:-1,:]+lon0[1:,:])/2.)   
    lat0 = np.squeeze((lat0[0,:,0:-1]+lat0[0,:,1:])/2.)   
    lat0 = np.squeeze((lat0[0:-1,:]+lat0[1:,:])/2.)   
    time0= np.squeeze(time0[:,0,0])

    return velDiv0,lon0,lat0,time0
###########################################################################################
def getCurl(period = [2,29],latlim=[23.0,31.0],lonlim=[-91.,-82.0],zlevel=0):

    uvel,vvel,time_v,lon_v,lat_v=getOsprein(period,latlim,lonlim,[zlevel,zlevel])
    lat0,time0,lon0 = np.meshgrid(lat_v,time_v,lon_v)
    print 'Computing velocity curl'
    x0 = (lon0-lon_v[0])*(40000./360.)*np.cos(np.deg2rad(lat0))
    y0 = (lat0-lat_v[0])*(40000./360.)

    dvdx=np.diff(np.squeeze(vvel),axis=2)/np.diff(x0,axis=2)
    dvdx = (dvdx[:,0:-1,:]+dvdx[:,1:,:])/2.
    dudy=np.diff(np.squeeze(uvel),axis=1)/np.diff(y0,axis=1)
    dudy = (dudy[:,:,0:-1]+dudy[:,:,1:])/2.
    velCurl0 = dvdx - dudy
    lon0 = np.squeeze((lon0[0,:,0:-1]+lon0[0,:,1:])/2.)   
    lon0 = np.squeeze((lon0[0:-1,:]+lon0[1:,:])/2.)   
    lat0 = np.squeeze((lat0[0,:,0:-1]+lat0[0,:,1:])/2.)   
    lat0 = np.squeeze((lat0[0:-1,:]+lat0[1:,:])/2.)   
    time0= np.squeeze(time0[:,0,0])

    return velCurl0,lon0,lat0,time0
###########################################################################################
def getStrain(period = [2,29],latlim=[23.0,31.0],lonlim=[-91.,-82.0],zlevel=0):

    uvel,vvel,time_v,lon_v,lat_v=getOsprein(period,latlim,lonlim,[zlevel,zlevel])
    lat0,time0,lon0 = np.meshgrid(lat_v,time_v,lon_v)
    print 'Computing strain'
    x0 = (lon0-lon_v[0])*(40000./360.)*np.cos(np.deg2rad(lat0))
    y0 = (lat0-lat_v[0])*(40000./360.)

    dvdx=np.diff(np.squeeze(vvel),axis=2)/np.diff(x0,axis=2)
    dvdx = (dvdx[:,0:-1,:]+dvdx[:,1:,:])/2.
    dudy=np.diff(np.squeeze(uvel),axis=1)/np.diff(y0,axis=1)
    dudy = (dudy[:,:,0:-1]+dudy[:,:,1:])/2.
    strain0 = dvdx + dudy
    lon0 = np.squeeze((lon0[0,:,0:-1]+lon0[0,:,1:])/2.)   
    lon0 = np.squeeze((lon0[0:-1,:]+lon0[1:,:])/2.)   
    lat0 = np.squeeze((lat0[0,:,0:-1]+lat0[0,:,1:])/2.)   
    lat0 = np.squeeze((lat0[0:-1,:]+lat0[1:,:])/2.)   
    time0= np.squeeze(time0[:,0,0])

    return strain0,lon0,lat0,time0
####################################################################################################
def OkuboWeiss(period = [2,29],latlim=[23.0,31.0],lonlim=[-91.,-82.0],zlevel=0):

    uvel,vvel,time_v,lon_v,lat_v=getOsprein(period,latlim,lonlim,[zlevel,zlevel])
    lat0,time0,lon0 = np.meshgrid(lat_v,time_v,lon_v)
    print 'Computing Okubo-Weiss parameter'
    x0 = (lon0-lon_v[0])*(40000./360.)*np.cos(np.deg2rad(lat0))
    y0 = (lat0-lat_v[0])*(40000./360.)

    dudx=np.diff(np.squeeze(uvel),axis=2)/np.diff(x0,axis=2)
    dudx = (dudx[:,0:-1,:]+dudx[:,1:,:])/2.

    dvdx=np.diff(np.squeeze(vvel),axis=2)/np.diff(x0,axis=2)
    dvdx = (dvdx[:,0:-1,:]+dvdx[:,1:,:])/2.

    dudy=np.diff(np.squeeze(uvel),axis=1)/np.diff(y0,axis=1)
    dudy = (dudy[:,:,0:-1]+dudy[:,:,1:])/2.

    dvdy=np.diff(np.squeeze(vvel),axis=1)/np.diff(y0,axis=1)
    dvdy = (dvdy[:,:,0:-1]+dvdy[:,:,1:])/2.

    curl = dvdx - dudy
    Sn = dudx - dvdy
    Ss = dvdx + dudy
    S2 = Sn **2 + Ss**2
    OW = S2 - curl**2
    lon0 = np.squeeze((lon0[0,:,0:-1]+lon0[0,:,1:])/2.)   
    lon0 = np.squeeze((lon0[0:-1,:]+lon0[1:,:])/2.)   
    lat0 = np.squeeze((lat0[0,:,0:-1]+lat0[0,:,1:])/2.)   
    lat0 = np.squeeze((lat0[0:-1,:]+lat0[1:,:])/2.)   
    time0= np.squeeze(time0[:,0,0])

    return OW,lon0,lat0,time0



###########################################################################################
def getEke(period = [2,29],latlim=[23.0,31.0],lonlim=[-91.,-82.0],zlevel=0):

    uvel,vvel,time,lon,lat=getOsprein(period,latlim,lonlim,[zlevel,zlevel])
    print 'Computing eddy kinetic energy'
    eke = np.squeeze((uvel**2. + vvel**2.)/2.)
    lon,lat = np.meshgrid(lon,lat)
    return eke,lon,lat,time
#############################################################################################
def timeInterpField(time0,Var0,time1,method='linear'):
    if np.ndim(Var0)==1:
       Var0 = np.reshape(Var0,[Var0.size,1,1])
    NX = np.size(Var0,2) #lon0.size
    NY = np.size(Var0,1) #lat0.size
    NT = np.size(time1,0)
    Var1=np.zeros(((NT,NY,NX)))
    print NT,NY,NX
    print np.size(time0)
    print np.ndim(time0)
    # interpolate in time
    for j in range(NY):
        for i in range(NX):
            #print j,i
            Var1[:,j,i] = interpolate.griddata(time0,np.squeeze(Var0[:,j,i]),time1, method=method)
    return np.squeeze(Var1)    
#############################################################################################
def spaceInterpField(lat0,lon0,Var0,lat1,lon1,method='cubic'):
    # interpolate in space
    # lat1 and lon1 are matrices
    # [time,drifter] or [time] in the case of center of mass
    if np.ndim(lat0)> 1:
       lat0 = np.squeeze(lat0[:,0])
       lon0 = np.squeeze(lon0[0,:])

    if np.ndim(lat1)<=1:
       lat1 = np.reshape(lat1,(np.size(lat1),1))
       lon1 = np.reshape(lon1,(np.size(lon1),1))
    NT = np.size(lat1,0)
    Var1 = np.zeros((NT,np.size(lat1,1)))

    for t in range(NT): 
        #print str(t)+'/'+str(NT)
        f = interpolate.interp2d(lon0, lat0, np.squeeze(Var0[t,:,:]), kind=method)
        for i in range(lon1[t,:].size):
            Var1[t,i] = f(np.squeeze(lon1[t,i]), np.squeeze(lat1[t,i]))
    return np.squeeze(Var1)
#############################################################################################
def animate_Var(Var='Div',period = [2,30],latlim=[-999,-999],lonlim=[-999,-999]):

### be sure of the starting day of the simulation
### in this case ([2,12]), the simulation started on jan 02 2013.  
     
    if (Var=='Div')|(Var=='div'):

       command = 'velVar0,lon0,lat0,time = getDiv(period'
       Vmin = -0.08 #velDiv0.min() 
       Vmax = 0.08 #velDiv0.min()*(-1)
       fname = 'velDiv.mp4'
       cname = 'seismic'

    elif (Var=='Curl')|(Var=='curl'):
       command = 'velVar0,lon0,lat0,time = getCurl(period'
       Vmin = -0.08 #velDiv0.min() 
       Vmax = 0.08 #velDiv0.min()*(-1)
       fname = 'velCurl.mp4'
       cname = 'seismic'
  
    elif (Var=='Strain')|(Var=='strain'):
       command = 'velVar0,lon0,lat0,time = getStrain(period'
       Vmin = -0.08 #velDiv0.min() 
       Vmax = 0.08 #velDiv0.min()*(-1)
       fname = 'velStrain.mp4'
       cname = 'seismic'

    elif (Var=='Eke')|(Var=='eke'):
       command = 'velVar0,lon0,lat0,time = getEke(period'
       Vmin = 0. #velDiv0.min() 
       Vmax = 0.1 #velDiv0.min()*(-1)
       fname = 'velEke.mp4'
       cname = 'jet'

    elif (Var=='Salinity')|(Var=='salinity'):
       command = 'velVar0,lon0,lat0,time = getSalinity(period'
       Vmin = 33 #velDiv0.min() 
       Vmax = 36 #velDiv0.min()*(-1)
       fname = 'velSal.mp4'
       cname = 'gist_rainbow'

    elif (Var=='Okubo-Weiss'): #|(Var=='strain'):
       command = 'velVar0,lon0,lat0,time = OkuboWeiss(period'
       Vmin = -0.05 #velDiv0.min() 
       Vmax = 0.05 #velDiv0.min()*(-1)
       fname = 'OkuboWeiss.mp4'
       cname = 'seismic'


# Plot limits
    global lamin, latmax, lonmin, lonmax, lon, lat, velVar
    if latlim[0] == -999:
       command = command+')'
    else:
       command = command+',latlim,lonlim)'

    exec command in globals(), locals()

    lon=lon0
    lat = lat0
 #   if (Var=='Salinity')|(Var=='salinity'):
 #      lon,lat = np.meshgrid(lon,lat) 
    velVar=velVar0  
    latmin = np.min(lat) 
    latmax = np.max(lat)
    lonmin = np.min(lon) 
    lonmax = np.max(lon)  
       
## gridcolor
    N = 200. #number of contour lines
    dL   = (Vmax-Vmin)/N
    V = np.arange(Vmin,Vmax+dL,dL) #colorscale
    figW = 20.
    figH = 20.
    FS = 42
    fig = pl.figure(figsize=(figW,figH))
    plot=fig.add_subplot(111)
    m = Basemap(llcrnrlat=latmin, urcrnrlat=latmax, llcrnrlon=lonmin,urcrnrlon=lonmax,  projection='cyl', resolution='h')   
    m.drawcoastlines()
    m.fillcontinents('0.7')
    parallels = np.arange(np.ceil(latmin),np.ceil(latmax),.5)
    m.drawparallels(parallels,labels=[True,False,False,False],fontsize=FS)
    meridians = np.arange(np.ceil(lonmin),np.ceil(lonmax),.5)
    m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=FS)
    date_text = plot.text(0.02, 1.02, '', transform=plot.transAxes,fontsize=FS,fontweight='bold')

    nframes = int(np.size(time))
    Day = np.floor(time-time[0])+1   
    Hour = np.floor(np.mod(time-time[0]+1,Day)*24)         

    def animate(j):
        im = m.contourf(lon,lat,np.squeeze(velVar[j,:,:]),V,cmap=plt.get_cmap(cname))
        cbar = m.colorbar(im,location='right',pad="5%")
#        cbar.set_label(cbarText,fontsize = FS)
        cbar.ax.tick_params(labelsize = FS)
        print nframes,j
        date_text.set_text('Day: ' + str(int(Day[j]))+', hour: '+ str(int(Hour[j])))
        return im,date_text
        
#    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=nframes, interval=5, blit=True)
    anim = animation.FuncAnimation(fig, animate, frames=nframes, interval=5, blit=True)
    fps=10
    frame_prefix='_tmp'
    anim.save(fname, fps=fps, codec='mpeg4')
    pl.show()
    return anim






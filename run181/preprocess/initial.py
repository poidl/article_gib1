# -*- coding: utf-8 -*-
#from Scientific.IO.NetCDF import NetCDFFile as NF
from netCDF4 import Dataset as NF
import matplotlib.pyplot as plt
import numpy as np
import mod_data as md
import os as os
import inspect as inspect

rundir= os.path.dirname(inspect.getfile(inspect.currentframe()))
grdname=rundir+'/../data/grd.nc';
initname=rundir+'/../data/init.nc';
bdyname=rundir+'/../data/bdy.nc';


grid=NF(grdname,'r')
h=grid.variables['D'][:]
lon=grid.variables['xh'][:]
lat=grid.variables['yh'][:]
grid.close()

#fac=np.linspace(0.,-0.5,h.shape[1])
#eta1=-np.ones(h.shape)*np.tile(fac,(h.shape[0],1)) 
eta1=np.zeros(h.shape)
e2east=80.
eta2a=eta1+e2east; eta2b=h-1e-9; #eta2b[eta2b<eta2a]=eta2a[eta2b<eta2a]+1e-9  
fac=np.linspace(1,0,10); mat=np.ones(h.shape); 
[x1,x2]=np.meshgrid(fac,np.arange(mat.shape[0]));
#mat[:,22:34]= x1;
#mat[:,34:]=0.;
mat[:,44:54]= x1;
mat[:,54:]=0.;
eta2=eta2b*mat+eta2a*(1-mat);
#eta2=eta1+10.
### rest
#eta2=75*np.ones(eta2.shape)
###
eta3=h;
eta=np.array([eta1,eta2,eta3])
eta_ini=-eta;

u=np.squeeze(0.*eta[0:2,...]);

lay=[1027.6,1030.1]

md.create_init(initname,grdname,'x')
init=NF(initname,'a')
init.variables['LAYER'][:]=lay;
init.variables['ETA'][:]=eta_ini
init.variables['u'][:]=u
init.variables['v'][:]=u
init.close()
eta_bdy=-eta;

md.create_bdy(bdyname,grdname,'x')
bdy=NF(bdyname,'a')
bdy.variables['LAYER'][:]=lay;
bdy.variables['ETA'][:]=eta_ini
bdy.variables['u'][:]=u
bdy.variables['v'][:]=u
bdy.close()

plt.close()
plt.close()
plt.close()
plt.close()
plt.close()
plt.close()
plt.close()
plt.figure(1)
ax = plt.axes()
#p=ax.contourf(lon,lat,eta_ini[1,...],20)
p=ax.contourf(eta_ini[1,...],20)
plt.colorbar(p)
#plt.xlim(lon[0], lon[-1])
#plt.ylim(lat[0], lat[-1])



plt.figure(2)
ax = plt.axes()
p=ax.imshow(mat,interpolation='nearest')
plt.colorbar(p)

plt.figure(3)
plt.plot(eta_ini[1,8,])

plt.figure(4)
plt.imshow(h,interpolation='nearest')
plt.colorbar()
plt.figure(5)
plt.imshow(eta_ini[2,...],interpolation='nearest')
plt.colorbar()

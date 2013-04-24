#from Scientific.IO.NetCDF import NetCDFFile as NF
from netCDF4 import Dataset as NF
import matplotlib.pyplot as plt
import numpy as np
import mod_data as md
import os as os
import inspect as inspect

rundir= os.path.dirname(inspect.getfile(inspect.currentframe()))
grdname=rundir+'/../data/grd.nc';
spngname=rundir+'/../data/sponge.nc';

grid=NF(grdname,'r')
h=grid.variables['D'][:]
grid.close()

idamp=np.zeros(h.shape)

#dmp=1./77.
dmp=1e-4
#dmp=0.;
#dmp=1./20;
n=15;
nv=dmp*np.linspace(1,0,n)
[x1,x2]=np.meshgrid(nv,np.arange(h.shape[0]))
idamp[:,:n]=x1

dmp=2.5e-5
n=150;
nv=dmp*np.linspace(1,0,n)
[x1,x2]=np.meshgrid(nv,np.arange(h.shape[0]))
idamp[:,-n:]=x1[:,-1::-1]
#idamp[n:-n,:n]=x1[n:-n]

#[sq,tr]=np.meshgrid(nv,np.arange(n))
#co=0
#for i in np.arange(n):
#    for j in np.arange(1+co,n):
#        sq[i,j]=sq[i,co]
#    co=co+1    
#
#idamp[:n,:n]=sq  
#idamp[-n:,:n]=sq[-1::-1,:]  
#
#x2=x1.transpose()
#ul=5 #upper left
#idamp[:n,n:ul+n]=x2[:,:ul]
#ll=12 #lower left
#idamp[-n:,n:ll+n]=x2[-1::-1,:ll]
  
md.create_spng(spngname,grdname,'x')
spng=NF(spngname,'a')
spng.variables['Idamp'][:]=idamp;
spng.close() 
  
#mask = np.load('mask.npy')
#maskp=np.ones(np.array(mask.shape))
#maskp[mask]=1.
#maskp[~mask]=0.

plt.close()
plt.close()
plt.close()
plt.close()
plt.close()
plt.figure(1)
ax = plt.axes()
p=ax.imshow(idamp,interpolation='nearest')
#c=ax.contour(maskp[-1::-1,:])
plt.colorbar(p)

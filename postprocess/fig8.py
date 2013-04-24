# -*- coding: utf-8 -*-
from netCDF4 import Dataset as NF
import numpy as np
import matplotlib.pyplot as plt


def volume_flux(sname,gname,tind,lonind):

	ff=NF(sname,'r')
	u=ff.variables['u'][-tind:,:,:,lonind:lonind+2]
	thick=ff.variables['h'][-tind:,:,:,lonind]
	time=ff.variables['Time'][-tind:]
	ff.close()

	ff=NF(gname,'r')
	dyh=ff.variables['dyh'][:,lonind].flatten()
	ff.close()
	u=u[:,:,:,:-1]+0.5*np.diff(u,axis=3) #regrid
	u=u[...,0]
	thick=thick[...,0]

	dyh=np.tile(dyh,(u.shape[0],u.shape[1],1))

	U=thick*u*dyh
	U=np.sum(U,axis=2)

	U=U[:,:]/1e6
	U1m=np.mean(U[:,0])
	U2m=np.mean(U[:,1])
	U1std=np.std(U[:,0])
	U2std=np.std(U[:,1])

	return U1m, U2m, U1std, U2std

basic=range(148,152)
dense=range(178,182)
light=range(182,186)

Dn=880.
tind=120

drho=[1.5,2.,2.5]
#for ii in range(4):
	#print str(ii) + ' of 3'
	#q=[]
	#qstd=[]
	#eta_exit=[]
	
	#for dr in [light,basic,dense]:

handles=[]
col=['b','g','r']
for dr,icol in zip([light,basic,dense],range(3)):
	print dr
	q=[]
	qstd=[]
	eta_exit=[]
	for ii in range(4):
	
		ddir='../run'+str(dr[ii])+'/';
		sname=ddir+'saves/save0.00e00.517.085.nc';
		aname=ddir+'saves/avfld0.00e00.517.085.nc';
		dname=ddir+'saves/D.517.85.2.nc';
		gname=ddir+'saves/grid.517.85.nc';
		
		ff=NF(dname,'r')
		D=ff.variables['D'][:]
		ff.close()
		
		jmid=np.argmax(D[:,0])
		inarrow=np.argwhere(D[jmid,:]==Dn)[0]
		iexit=np.argwhere(D[10,inarrow:]>1)[0]+inarrow-1 #index of strait exit
		
		
		
		ff=NF(aname,'r')
		etm=ff.variables['etm'][-tind:,1,jmid,:]
		ff.close()
		
		em_exit=np.mean(etm[:,iexit],axis=0)
		U1m, U2m, U1std, U2std = volume_flux(sname,gname,tind,inarrow)
		q.append(U1m)
		qstd.append(U1std)
		eta_exit.append(em_exit)
		
	plt.plot(eta_exit,q,'*',markersize='10',color=col[icol])
	p=plt.plot(eta_exit,q,color=col[icol])
	handles.append(p[0])
	
leg = plt.legend(handles,('LIGHT','BASIC','DENSE'),loc=(0.7,0.06))
plt.xlabel('$\eta_{exit}$'+' at y=0 km')
plt.ylabel('$Q_1$')
plt.grid('on')    
plt.savefig('figures/fig8.pdf')


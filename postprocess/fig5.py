# -*- coding: utf-8 -*-
from netCDF4 import Dataset as NF
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

plt.close('all')

tind=120
ifc=1
lay=1
Ds=284.
Dn=880.
runnr=150

vec=range(148,152)

#vec=range(166,170)

ddir='../run148/'
gname=ddir+'saves/grid.517.85.nc';
gname1=ddir+'data/grd.nc';
ininame=ddir+'data/init.nc';
dname=ddir+'saves/D.517.85.2.nc';

ff=NF(ininame,'r')
xh=ff.variables['xh'][:]*1e-3
yh=ff.variables['yh'][:]*1e-3
ff.close()

ff=NF(gname1,'r')
D=ff.variables['D'][:]
ff.close()

ff=NF(dname,'r')
r1=ff.variables['R'][0]
r2=ff.variables['R'][1]
ff.close()

jmid=np.argmax(D[:,0])
isill=np.argwhere(D[jmid,:]==Ds)[0]
inarrow=np.argwhere(D[jmid,:]==Dn)[0]
iexit=np.argwhere(D[10,inarrow:]>1)[0]+inarrow-1 #index of strait exit

ff=NF(gname,'r')
wet=ff.variables['wet'][:]
dyh=ff.variables['dyh'][:]
ff.close()

left=isill-2; right=isill+38


fac=0.8
fig=plt.figure(figsize=(fac*13,fac*8))
col=['r','b','g','c']
handles=[]

for jj in range(len(vec)):
#for jj in [0]:
	print 'processing '+str(jj)+' of '+str(len(vec))
	ddir='../run'+str(vec[jj])+'/';
	sname=ddir+'saves/save0.00e00.517.085.nc';
	aname=ddir+'saves/avfld0.00e00.517.085.nc';

	ff=NF(sname,'r')
	thick=ff.variables['h'][-tind:,:2,:,:]
	u=ff.variables['u'][-tind:,:2,:,:]
	ff.close()

	ff=NF(aname,'r')
	e=ff.variables['etm'][-tind:,ifc,:,:]
	e0=ff.variables['etm'][-tind:,0,:,:]
	ff.close()

	u=u[:,:,:,:-1]+0.5*np.diff(u,axis=3) #regrid
	gd=9.81*(r2-r1)/r2
	fr=u**2/(gd*thick)
	fr_comp=fr[:,0,:,:]+fr[:,1,:,:]
	fr_comp=np.mean(fr_comp,axis=0)
	fr_upper=np.mean(fr[:,0,:,:],axis=0)
	em=np.mean(e,axis=0)
	
	
	####### Fr triangular cross section
	## Bormans (nonrectangular)
	#eta_1=e0[:,jmid,:]
	#W=np.sum( wet*dyh, axis=0 ); W=np.array([W]*eta_1.shape[0])
	#b=-D[jmid,:]; b=np.array([b]*eta_1.shape[0])
	#H=eta_1-b
	#alpha=W/(2*H)
	#eta_2=e[:,jmid,:]
	#A_1=W*(eta_1-eta_2)-alpha*(eta_1-eta_2)**2.
	#A_2=alpha*(eta_2-b)**2.
	#W_int=2.*alpha*(eta_2-b)
	u_1=u[:,0,jmid,:]
	u_2=u[:,1,jmid,:]
	#u_1=np.mean( u[:,0,jmid-2:jmid+3,:],axis=1 )
	#u_2=np.mean( u[:,1,jmid-2:jmid+3,:],axis=1 )
	
	#f_borm=u_1**2*W_int/(gd*A_1) + u_2**2*W_int/(gd*A_2)
	#f_borm=np.mean(f_borm,axis=0)
	
	# Henderson 66 (Bryden and Kinder 91, p.453)
	# identical to f_borm but simpler expression
	H1s=-e[:,jmid,:]
	D_bk=D[jmid,:]; D_bk=np.array([D_bk]*H1s.shape[0])
	h1s=H1s*(D_bk-0.5*H1s)/(D_bk-H1s)
	h2s=0.5*(D_bk-H1s)
	
	f_bk=u_1**2/(gd*h1s)+u_2**2/(gd*h2s)
	f_bk=np.mean(f_bk,axis=0)
	

	#clip 
	fr_comp_mean=np.mean(fr_comp[jmid-2:jmid+3,left:right],axis=0)
	fr_comp=fr_comp[jmid,left:right]
	fr_upper=fr_upper[jmid,left:right]
	#f_borm=f_borm[left:right]
	f_bk=f_bk[left:iexit+1]
	
	em_mean=np.mean(em[jmid-2:jmid+3,left:right],axis=0)
	em=em[jmid,left:right]

	a2=[0.08,0.07,0.9,0.4]
	a1=[0.08,0.55,0.9,0.4]
	ax1 = plt.axes(a1)
	p=plt.plot(xh[left:right],fr_comp,col[jj])
	
	handles.append(p[0])
	#plt.plot(xh[left:right],fr_comp_mean,col[jj]+'--')
	#plt.plot(xh[left:iexit+1],f_bk,col[jj]+'_')

	plt.plot(xh[left:iexit+1],f_bk, marker='v', linestyle='None',  color=col[jj], markersize=6)
	plt.plot(xh[left:right],fr_upper,marker='o',linestyle='None',color=col[jj],markersize=8,markerfacecolor='none',markeredgecolor=col[jj])
	
	ax2=plt.axes(a2)
	plt.plot(xh[left:right],em,col[jj])
	#plt.plot(xh[left:right],em_mean,col[jj]+'--')	



ax1=plt.axes(a1)		
for ii in [isill,inarrow,iexit]:
	plt.axvline(xh[ii],ymin=0,ymax=1,color='k',linewidth=2)
plt.axhline(1,xmin=0,xmax=1,color='k',linewidth=2)
plt.xlim([-12,65])
plt.ylim([0,3.9])
#plt.xlabel('x [km]',axes=ax1)
plt.ylabel('$G^2,$'+' '+'$G_{triang}^2,$'+' '+'$F_1^2$',fontsize=16)
plt.title('Froude Number')
leg = ax1.legend(handles,('BASIC-20','BASIC-40','BASIC-60','BASIC-80'),loc=(0.1,0.4))
plt.text(-0.08,1,'a)',transform = ax1.transAxes,fontsize=15)
plt.grid('on')

ax2=plt.axes(a2)
#plt.plot(xh[left:right],-bot,color='k')
for ii in [isill,inarrow,iexit]:
	plt.axvline(xh[ii],ymin=0,ymax=1,color='k',linewidth=2)
plt.xlim([-12,65])
plt.ylim([-155,-35])
plt.xlabel('x [km]')
plt.ylabel('Depth [m]')
plt.title('Interface Depth $\eta$ [m]')
plt.text(-0.08,1,'b)',transform = ax2.transAxes,fontsize=15)
plt.grid('on')
for s,ii in zip(['Sill','Contraction','Exit'],[isill,inarrow,iexit]):
	plt.annotate(s, [xh[ii],-140.], xytext=[xh[ii],-175.], xycoords='data',
	 arrowprops=dict(arrowstyle="-",facecolor='black'),textcoords='data',horizontalalignment='center')


plt.subplots_adjust(left=0.0, bottom=0.01, right=1.0, top=0.98, wspace=0.0, hspace=0.05)

plt.savefig('figures/fig5.pdf')






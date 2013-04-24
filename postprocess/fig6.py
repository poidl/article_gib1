# -*- coding: utf-8 -*-
from netCDF4 import Dataset as NF
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from matplotlib.path import Path
import matplotlib.patches as patches

plt.close('all')

tind=120
ifc=1
lay=0
Ds=284.
Dn=880.
runnr=150

vec=range(148,152)

#vec=range(166,170)

fac=1.
fig=plt.figure(figsize=(fac*13,fac*7))

g=[]
cblevs=[]

for jj in range(len(vec)):
        print 'processing '+str(jj)+' of '+str(len(vec))
	ddir='../run'+str(vec[jj])+'/';
	sname=ddir+'saves/save0.00e00.517.085.nc';
	aname=ddir+'saves/avfld0.00e00.517.085.nc';
	gname=ddir+'saves/grid.517.85.nc';
	gname1=ddir+'data/grd.nc';
	ininame=ddir+'data/init.nc';
	dname=ddir+'saves/D.517.85.2.nc';

	ff=NF(ininame,'r')
	xh=ff.variables['xh'][:]*1e-3
	yh=ff.variables['yh'][:]*1e-3
	ff.close()

	ff=NF(sname,'r')
	thick=ff.variables['h'][-tind:,lay,:,:]
	u=ff.variables['u'][-tind:,lay,:,:]
	ff.close()

	ff=NF(aname,'r')
	e=ff.variables['etm'][-tind:,ifc,:,:]
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
	ff.close()

	u=u[:,:,:-1]+0.5*np.diff(u,axis=2) #regrid
	gd=9.81*(r2-r1)/r2
	#fr=u**2/(gd*thick) (armi & farmer)
	fr=u/np.sqrt(gd*thick)
	fr=np.mean(fr,axis=0)
	em=np.mean(e,axis=0)
        #em[D+em<1e-2]=np.nan
        bo=(D+em<1e-1)

	#clip 
	left=isill-3; right=isill+38
        bottom=jmid-8; top=jmid+9
	xh=xh[left:right]
        yh=yh[bottom:top]
        fr=fr[bottom:top,left:right]
	em=em[bottom:top,left:right]
	wet=wet[bottom:top,left:right]
        bo=bo[bottom:top,left:right]
	#plt.figure(figsize=(25,10))

	ax=fig.add_subplot(2,2,jj+1, aspect='equal')
	#cmax=np.nanmax(np.abs(em))
	#cmin=-cmax
        if ifc==1:
		cmin=0.
		cmax=np.nanmax(fr)
		if cmax>4.:
		  cmax=4.
	fr[wet==0]=np.nan
	palette = matplotlib.cm.jet
	norm1=matplotlib.colors.Normalize(vmin = cmin, vmax = cmax, clip = False)
        do=1./np.arange(15)[1:] # hit 1
        do2=(cmax-cmin)/14
        dl=do[np.argmin(np.abs(do-do2))]
        levs=np.arange(cmin,cmax,dl)
        #dlev=np.diff(levs)[0]
        if np.argwhere(levs[::2]==1.):
	  cblevs.append( levs[::2] )
	else:
	  cblevs.append( levs[1::2] )
	gt=plt.contourf(xh,yh,fr,levs,\
	    cmap=palette,\
	    norm = norm1, extend='both'
	    )
	    
	g.append(gt)
	
	con=plt.contour(xh,yh,fr,levels=[1.],colors='k',linestyles='-')
        for kk in range(bo.shape[0]):
		for ll in range(bo.shape[1]):
			if bo[kk,ll]==True:
        			plt.plot(xh[ll],yh[kk],marker='x',color='white',markersize=8)

######################
	dx2=np.diff(xh)[0]*0.5

	iy=[]
	for ii in range(wet.shape[1]):
		li=np.where(wet[:,ii]==1)[0][0]
		iy.append(li)
		
	verts = [
		(xh[0], yh[0]), # left, bottom
		(xh[-1], yh[0]), # right, bottom
		]
	codes = [Path.MOVETO,Path.LINETO]

	for ii in range(wet.shape[1]-1):
		if iy[ii]<iy[ii+1]:
			iy[ii]=iy[ii+1]
	
	for ii in range(wet.shape[1])[-1::-1]:
		if ii<wet.shape[1]-1:
			if iy[ii]!=iy[ii+1]:
				verts.append( (xh[ii+1],yh[iy[ii]]) )
				codes.append( Path.LINETO )

		verts.append( (xh[ii],yh[iy[ii]]) )
		codes.append( Path.LINETO )

	verts.append( (xh[0], yh[0]) ) # ignored
	codes.append( Path.CLOSEPOLY )

	path = Path(verts,codes)
	patch = patches.PathPatch(path, facecolor=(0.8,0.8,0.8), edgecolor='none')
	ax.add_patch(patch)
	
	nverts=[]
	for ii in verts:
		nverts.append( (ii[0],-1*ii[1]) )
	path = Path(nverts,codes)
	patch = patches.PathPatch(path, facecolor=(0.8,0.8,0.8), edgecolor='none')
	ax.add_patch(patch)
	
	plt.show()
    
#######################	

	plt.xlabel('x [km]')
	plt.ylabel('y [km]')
	plt.grid(True)

        txt=['a','b','c','d']
	plt.text(-0.15,1,str(txt[jj])+')',transform = ax.transAxes,fontsize=15)
	
plt.subplots_adjust(left=0.05, bottom=0.02, right=0.99, top=0.95, wspace=0.1, hspace=0.1)

myformat = plt.FormatStrFormatter('%2.1f')
for jj,gt,cbl in zip(range(4),g,cblevs):
	ax=fig.add_subplot(2,2,jj+1, aspect='equal')
	plt.colorbar(gt,orientation='horizontal',format=myformat,\
                     ticks=cbl,extend='both')
plt.savefig('figures/fig6.pdf')





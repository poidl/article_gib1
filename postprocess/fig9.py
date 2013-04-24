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
Ds=284.
Dn=880.
runnr=150

#vec=range(148,152)
#vec=range(166,170)
#vec=range(178,182)
#vec=range(182,186)
vec=[178,182,179,183]
#vec=[166,148,167,149]
#vec=range(152,155)
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

	ff=NF(aname,'r')
	e=ff.variables['etm'][-tind:,ifc,:,:]
	ff.close()


	ff=NF(gname1,'r')
	D=ff.variables['D'][:]
	ff.close()

	ff=NF(dname,'r')
	r2=ff.variables['R'][1]
	ff.close()

	jmid=np.argmax(D[:,0])
	isill=np.argwhere(D[jmid,:]==Ds)[0]
	inarrow=np.argwhere(D[jmid,:]==Dn)[0]
	iexit=np.argwhere(D[10,inarrow:]>1)[0]+inarrow-1 #index of strait exit

	ff=NF(gname,'r')
	wet=ff.variables['wet'][:]
	ff.close()

	em=np.mean(e,axis=0)
        #em[D+em<1e-2]=np.nan
        bo=(D+em<1e-1)
        cmid=em[jmid,inarrow]
	#clip 
	left=isill-3; right=isill+38
        bottom=jmid-8; top=jmid+9
	xh=xh[left:right]
        yh=yh[bottom:top]
	em=em[bottom:top,left:right]
	wet=wet[bottom:top,left:right]
        bo=bo[bottom:top,left:right]
	#plt.figure(figsize=(25,10))

	ax=fig.add_subplot(2,2,jj+1, aspect='equal')
	#cmax=np.nanmax(np.abs(em))
	#cmin=-cmax
        if ifc==1:
		cmin=cmid-60-5
		cmax=cmid+60+5
	em[wet==0]=np.nan
	palette = matplotlib.cm.RdBu
	norm1=matplotlib.colors.Normalize(vmin = cmin, vmax = cmax, clip = False)
        levs=np.linspace(cmin,cmax,14)
        dlev=np.diff(levs)[0]
        cblevs.append( np.linspace(levs[0]+0.5*dlev,levs[-1]-0.5*dlev,5) )
	gt=plt.contourf(xh,yh,em,levs,\
	    cmap=palette,\
	    norm = norm1, extend='both'
	    )
	g.append(gt)
	
	#con=plt.contour(xh,yh,em,cblevs[1:-1],colors='k',linestyles='-')
        for kk in range(bo.shape[0]):
		for ll in range(bo.shape[1]):
			if bo[kk,ll]==True:
        			plt.plot(xh[ll],yh[kk],marker='x',color='green',markersize=8)

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
	#plt.title(r'$\eta$'+' [m]')
	plt.grid(True)
    
       

        txt=['a','b','c','d']
	plt.text(-0.15,1,str(txt[jj])+')',transform = ax.transAxes,fontsize=15)

plt.subplots_adjust(left=0.05, bottom=0.02, right=0.99, top=0.95, wspace=0.1, hspace=0.1)

myformat = plt.FormatStrFormatter('%2.1f')
for jj,gt,cbl in zip(range(4),g,cblevs):
	ax=fig.add_subplot(2,2,jj+1, aspect='equal')
	plt.colorbar(gt,orientation='horizontal',format=myformat,\
                     ticks=cbl,extend='both')
    
    
                     
plt.savefig('figures/fig9.pdf')





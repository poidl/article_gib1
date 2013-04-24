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
vec=range(148,152)
#vec=range(166,170)
#vec=range(178,182)
######vec=range(170,174) #cape
#vec=range(188,192) #fplane
vec=range(178,182) #dense
vec=range(182,186) #light
#vec=range(192,196) #cape-light
#vec=range(196,200) #cape-dense
#vec=range(200,204) #cape-fplane
#vec=[206]
#vec=range(207,211) #western basin 300m
#vec=range(211,215) # three layers
#vec=range(215,219) # three layers cape
#vec=range(219,223) #CAPE THREE LAYERS MW LARGEGRAD: like 215-218 but  [1027.6,1029.6,1029.8]
#vec=range(228,232)
#vec=range(232,236)
#####vec=range(236,240)
#vec=range(224,186)
vec=[237]
fac=1.0
fig=plt.figure(figsize=(fac*9.2,fac*3.2))

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
	bo=(D+em<1e-1)
        cmid=em[jmid,inarrow]
        
        embasin=em.copy()
	embasin[:,:iexit]=0.
	#clip sponge
	left=isill-3; right=150
	xh=xh[left:-right]
	em=em[:,left:-right]
	wet=wet[:,left:-right]
        bo=bo[:,left:-right]
        embasin=embasin[:,left:-right]
	#plt.figure(figsize=(25,10))
	

	jmin,imin=np.nonzero(embasin.min()==embasin)
        if jmin.shape[0]>1:
	  print 'warning: multiple minima'

	ax=fig.add_subplot(len(vec),1,jj+1, aspect='equal')
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
	con=plt.contour(xh,yh,em,levels=[-120],colors='k',linestyles='-')
	con=plt.contour(xh,yh,em,levels=[-130],colors='k',linestyles='--')
	con=plt.contour(xh,yh,em,levels=[-140],colors='k',linestyles=':')
	if 0:
		for kk in range(bo.shape[0]):
		  for ll in range(bo.shape[1]):
			  if bo[kk,ll]==True:
				  plt.plot(xh[ll],yh[kk],marker='x',color='green',markersize=1)
	print 'min(e): ' +str(em[jmin,imin])
	print 'xpos: ' +str(xh[imin])
	print 'ypos: ' +str(yh[jmin])
	#plt.plot(xh[imin],yh[jmin],marker='o',color='yellow',markersize=5)

######################
# south shore
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
	patch = patches.PathPatch(path, facecolor=(0.7,0.7,0.7), edgecolor='none')
	ax.add_patch(patch)
	
###
# north shore
	iy=[]
	for ii in range(wet.shape[1]):
		li=wet.shape[1]-np.where(wet[:,ii]==1)[0][-1]
		iy.append(li)
		
	verts = [
		(xh[0], yh[-1]), # left, top
		(xh[-1], yh[-1]), # right, top
		]
	codes = [Path.MOVETO,Path.LINETO]

	for ii in range(wet.shape[1]-1):
		if iy[ii]<iy[ii+1]:
			iy[ii]=iy[ii+1]
	
	for ii in range(wet.shape[1])[-1::-1]:
		if ii<wet.shape[1]-1:
			if iy[ii]!=iy[ii+1]:
				verts.append( (xh[ii+1],yh[wet.shape[1]- iy[ii]]) )
				codes.append( Path.LINETO )

		verts.append( (xh[ii],yh[wet.shape[1]-iy[ii]]) )
		codes.append( Path.LINETO )

	verts.append( (xh[0], yh[-1]) ) # ignored
	codes.append( Path.CLOSEPOLY )

	path = Path(verts,codes)
	patch = patches.PathPatch(path, facecolor=(0.7,0.7,0.7), edgecolor='none')
	ax.add_patch(patch)
	
	plt.show()
    
#######################

	plt.xlabel('x [km]')
	plt.ylabel('y [km]')
       # plt.title(r'$\overline{\eta}$'+' [m]')
                     
                     
        plt.grid(True)             
        txt=['a','b','c','d']
	#plt.text(-0.08,1,str(txt[jj])+')',transform = ax.transAxes,fontsize=15)

plt.subplots_adjust(left=0.0, bottom=0.01, right=1.06, top=0.98, wspace=0.0, hspace=0.05)

myformat = plt.FormatStrFormatter('%2.1f')
for jj,gt,cbl in zip(range(4),g,cblevs):
	ax=fig.add_subplot(len(vec),1,jj+1, aspect='equal')
	cb=plt.colorbar(gt,orientation='horizontal',format=myformat,\
			ticks=cbl,extend='both')
	#cb=plt.colorbar(gt,format=myformat,\
	#		ticks=cbl,extend='both')
                     
	bb=cb.ax.get_position()
	#cb.ax.set_position([pos[0][0], pos[0][1],pos[1][0],0.5*pos[1][1]])
	cb.ax.set_position([ bb.x0, bb.y0+0.05, bb.width, 0.5*bb.height])
plt.savefig('figures/fig12.pdf')





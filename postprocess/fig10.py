# -*- coding: utf-8 -*-
from netCDF4 import Dataset as NF
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from matplotlib.path import Path
import matplotlib.patches as patches

plt.close('all')


ifc=1
Ds=284.
Dn=880.
runnr=170
vec=range(1,7)
vec=range(1,37,7)
vec=range(2,30,5)
vec[3]=38
vec[4]=41
vec[5]=44
vec=[vec[0],vec[3],vec[1],vec[4],vec[2],vec[5]]
print vec
#vec=range(2,33,6)

#vec=range(182,186)
fac=1.0
fig=plt.figure(figsize=(fac*16,fac*10))
ax = fig.add_subplot(111,frame_on=False)
#ax.xaxis.set_ticks_position("none")
#ax.yaxis.set_ticks_position("none")
ax.set_xticks([]) 
ax.set_yticks([]) 

ddir='../run'+str(runnr)+'/';
sname=ddir+'saves/save0.00e00.517.085.nc';
aname=ddir+'saves/avfld0.00e00.517.085.nc';
gname=ddir+'saves/grid.517.85.nc';
gname1=ddir+'data/grd.nc';
ininame=ddir+'data/init.nc';
dname=ddir+'saves/D.517.85.2.nc';



for jj in range(len(vec)):
        print 'processing '+str(jj)+' of '+str(len(vec))
	ff=NF(ininame,'r')
	xh=ff.variables['xh'][:]*1e-3
	yh=ff.variables['yh'][:]*1e-3
	ff.close()

	ff=NF(gname1,'r')
	D=ff.variables['D'][:]
	ff.close()

	ff=NF(dname,'r')
	r2=ff.variables['R'][1]
	ff.close()

	ff=NF(gname,'r')
	wet=ff.variables['wet'][:]
	ff.close()
	
	ff=NF(sname,'r')
	e=ff.variables['e'][vec[jj],ifc,:,:]
	time=ff.variables['Time'][vec[jj]]
	ff.close()
	
        ff=NF(aname,'r')
	etm=ff.variables['etm'][vec[jj],ifc,:,:]
	ff.close()
	
	jmid=np.argmax(D[:,0])
	isill=np.argwhere(D[jmid,:]==Ds)[0]
	inarrow=np.argwhere(D[jmid,:]==Dn)[0]
	iexit=np.argwhere(D[10,inarrow:]>1)[0]+inarrow-1 #index of strait exit
	
	bo=(D+etm<1e-1)
        cmid=e[jmid,inarrow] 

	#clip sponge
	left=isill-3; right=250
	xh=xh[left:-right]
	e=e[:,left:-right]
	wet=wet[:,left:-right]
        bo=bo[:,left:-right]
        
	#plt.figure(figsize=(25,10))

	axs=fig.add_subplot(3,2,jj+1, aspect='equal')
	if ifc==1:
	  cmin=-140
	  cmax=-10

	e[wet==0]=np.nan
	bo[wet==0]=False
	palette = matplotlib.cm.RdBu
	norm1=matplotlib.colors.Normalize(vmin = cmin, vmax = cmax, clip = False)
        levs=np.linspace(cmin,cmax,14)
        #dlev=np.diff(levs)[0]
        #cblevs=np.linspace(levs[0]+0.5*dlev,levs[-1]-0.5*dlev,5)
        cblevs=levs
	g=plt.contourf(xh,yh,e,levs,\
	    cmap=palette,\
	    norm = norm1, extend='both'
	    )
	#con=plt.contour(xh,yh,e,levels=[-120],colors='k',linestyles='-')
	#con=plt.contour(xh,yh,e,levels=[-130],colors='k',linestyles='--')
	#con=plt.contour(xh,yh,e,levels=[-140],colors='k',linestyles=':')

        if 0: # horz. velocity arrows
	  ff=NF(sname,'r')
	  u=ff.variables['u'][vec[jj],0,...]
	  v=ff.variables['v'][vec[jj],0,...]
	  ff.close()
	  u=u[:,:-1]+0.5*np.diff(u,axis=1)
          v=v[:-1,:]+0.5*np.diff(v,axis=0)
          u=u[:,left:-right]
          v=v[:,left:-right]
          yg=np.arange(u.shape[0])
	  xg=np.arange(u.shape[1])
          isub=2
	  jsub=2
	  utmp=u[::jsub,::isub]
	  vtmp=v[::jsub,::isub]
	  xg=xg[::isub]
	  yg=yg[::jsub]
	  print xh[::isub].shape
	  print yh[::jsub].shape
	  print utmp.shape
	  qu=plt.quiver(xh[::isub],yh[::jsub],utmp,vtmp,scale=30.)
	  qk = plt.quiverkey(qu, 0.8, 1.03, 1, r'$1 m/s$', labelpos='E',
               fontproperties={'size':15,'weight': 'bold'})

	if 0: 
	  for kk in range(bo.shape[0]):
	    for ll in range(bo.shape[1]):
		    if bo[kk,ll]==True:
			    plt.plot(xh[ll],yh[kk],marker='x',color='green',markersize=1)

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
	axs.add_patch(patch)
	
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
	axs.add_patch(patch)
	
	plt.show()
    
#######################
	if jj in [4,5]:
		plt.xlabel('x [km]')
	if jj in [0,2,4]:
		plt.ylabel('y [km]')
        plt.title(r'$\eta$'+' [m], t='+str(time)+' d')

        myformat = plt.FormatStrFormatter('%2.1f')
        #plt.colorbar(g,format=myformat,\
        #             ticks=cblevs,extend='both')
                     
                     
        plt.grid(True)             
        txt=['a','d','b','e','c','f']
	if jj in [1,3,5]:
		plt.text(-0.04,1,str(txt[jj])+')',transform = axs.transAxes,fontsize=15)
	else:
		plt.text(-0.09,1,str(txt[jj])+')',transform = axs.transAxes,fontsize=15)

plt.subplots_adjust(left=0.04, bottom=0.02, right=0.928, top=0.98, wspace=0.058, hspace=0.025)

cax=fig.add_axes([ax.get_position().x1+0.02,ax.get_position().y0,0.015,ax.get_position().y1-ax.get_position().y0])
plt.colorbar(g,cax,ticks=cblevs,format=myformat)


plt.savefig('figures/fig10.pdf')





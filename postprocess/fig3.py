# -*- coding: utf-8 -*-
#from Scientific.IO.NetCDF import NetCDFFile as NF
import matplotlib.pyplot as plt
from netCDF4 import Dataset as NF
import numpy as np

def pl(nr,nc,kk,yh,e,sd_e,D,isill,inarrow,iexit,jmid,e_ini,r2):

    ax=plt.subplot(nr,nc,kk)
    for ii,jj in zip(['r','b','g'],[isill,inarrow,iexit]):
            plt.plot(yh,e[:,jj],ii,lw=2)
            #plt.errorbar(yh,e[:,jj].flatten(), xerr=[0.]*len(yh),yerr=sd_e[:,jj].flatten(),\
            #            fmt=None,capsize=20,ecolor='k')    
	    plt.axhline(e[jmid,jj], xmin=0,xmax=1,linestyle='--',color=ii,linewidth=2)
            plt.plot(yh,-D[:,jj],ii)
	    inds=np.argwhere(~np.isnan(D[:,jj].flatten())).flatten()
	    y1=-D[:,jj][inds[0]]
	    plt.plot([yh[inds[0]],yh[inds[0]]],[y1,0.],ii)
	    plt.plot([yh[inds[-1]],yh[inds[-1]]],[y1,0.],ii)
    txt=['a','b','c','d']
    plt.text(-0.1,1,str(txt[kk-1])+')',transform = ax.transAxes)
    #plt.title(str(e_ini)+',  '+str(r2))
    plt.xlabel('y [km]')
    plt.ylabel('Depth [m]')
    
    plt.grid(True)
    plt.ylim([-210,0])
    #plt.xlim([39,47])
 
tind=120
ifc=1
Ds=284.
Dn=880.
vec=range(148,152)
#vec=range(166,170)
#vec=range(178,182)
#vec=range(182,186)
#vec=[148,206, 149,150]
#vec=range(170,174)
nr=2;nc=2

plt.close('all')
fac=1.0
fig=plt.figure(figsize=(fac*13,fac*6))


ii=1
for jj in range(len(vec)):
        print 'processing '+str(ii)+' of '+str(len(vec))
	ddir='../run'+str(vec[jj])+'/';
	sname=ddir+'saves/save0.00e00.517.085.nc';
        aname=ddir+'saves/avfld0.00e00.517.085.nc';
        gname=ddir+'saves/grid.517.85.nc';
	gname1=ddir+'data/grd.nc';
        ininame=ddir+'data/init.nc';
        dname=ddir+'saves/D.517.85.2.nc';
        print dname

        ff=NF(ininame,'r')
	e_ini=ff.variables['ETA'][1,10,-1]
        yh=ff.variables['yh'][:]*1e-3
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

        ff=NF(sname,'r')
        e=ff.variables['e'][-tind:,ifc,:,:]
        ff.close()

        ff=NF(aname,'r')
        etm=ff.variables['etm'][-tind:,ifc,:,:]
        ff.close()
        
        em=np.mean(etm,axis=0)
        sd_e=np.std(e,axis=0)

        em[wet==0]=np.nan
        sd_e[wet==0]=np.nan

        #e[D+e<1e-2]=np.nan
	D[wet==0]=np.nan   

        pl(nr,nc,ii,yh,em,sd_e,D,isill,inarrow,iexit,jmid,e_ini,r2)
        ii+=1

plt.savefig('figures/fig3.pdf')




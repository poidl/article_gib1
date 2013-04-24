# -*- coding: utf-8 -*-
from netCDF4 import Dataset as NF
import numpy as np


def volume_flux(sname,gname,tind,lonind):
        print lonind
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


tind=120
Ds=284.
Dn=880.
vec=range(148,152)
#vec=range(166,170)
#vec=range(170,174)
#vec=range(207,211)
#vec=range(178,182)
#vec=range(182,186)
#vec=range(215,219)
#vec=range(219,223)

stats=[[] for i in range(len(vec))]
for jj in range(len(vec)):
	ddir='../run'+str(vec[jj])+'/';
	sname=ddir+'saves/save0.00e00.517.085.nc';
	statsname=ddir+'saves/timestats.nc';
	aname=ddir+'saves/avfld0.00e00.517.085.nc';
        dname=ddir+'saves/D.517.85.2.nc';
        gname=ddir+'saves/grid.517.85.nc';
	grdname=ddir+'data/grd.nc';
	spngname=ddir+'data/sponge.nc';
        ininame=ddir+'data/init.nc';

        ff=NF(dname,'r')
        D=ff.variables['D'][:]
        ff.close()

        jmid=np.argmax(D[:,0])
        isill=np.argwhere(D[jmid,:]==Ds)[0]
        inarrow=np.argwhere(D[jmid,:]==Dn)[0]
        iexit=np.argwhere(D[10,inarrow:]>1)[0]+inarrow-1 #index of strait exit

        ff=NF(ininame,'r')
	e_ini=ff.variables['ETA'][1,10,-1]
        yh=ff.variables['yh'][:]*1e-3
	ff.close()

        ff=NF(sname,'r')
        e=ff.variables['e'][-tind:,1,jmid,:]
        ff.close()

        ff=NF(aname,'r')
        etm=ff.variables['etm'][-tind:,1,jmid,:]
        ff.close()

        
        em=np.mean(etm[:,[isill,inarrow,iexit]],axis=0).flatten()
        sd_e=np.std(e[:,[isill,inarrow,iexit]],axis=0).flatten()
        U1m, U2m, U1std, U2std = volume_flux(sname,gname,tind,inarrow)
        #print stats
	#print ' '
        stats[jj]=[em[0],em[1],em[2],U1m, U1m,\
			sd_e[0],sd_e[1],sd_e[2], U1std]

for st in stats:
	st[4]=(st[4]/stats[0][3]) #ratio maximal exchange

s=[]; ncol=10
s.append( ' \\begin{table}\\scriptsize ' )
s.append( ' \centering ' )
s.append( ' \\begin{tabular}{'+ 'c'*ncol +'} ' )
s.append( ' \\toprule ' )
s.append( ' &\multicolumn{4}{c}{Mean}& $\overline{Q_1}/\overline{Q_1}^{max}$ &\multicolumn{4}{c}{Standard Deviation}\\\\ ' )
s.append( ' \cmidrule{2-10} ' )
s.append( ' & \multicolumn{3}{c}{$\eta$} & $Q_1$ &  & \multicolumn{3}{c}{$\eta$} & $Q_1$ \\\\ ')
s.append( ' \cmidrule{2-10} ' )
s.append( ' &Sill&Contraction&Exit& & &Sill&Contraction&Exit&  \\\\ ' )
s.append( ' \cmidrule{2-10} ' )
label=['BASIC-20','BASIC-40','BASIC-60','BASIC-80']
def sstr(x,p):
	fs='%.'+str(int(p))+'f'
	return fs% x
precision=[1,1,1,3,3,  2,2,2,3]
for ii in range(len(label)):
	sl=label[ii]
	for st,jj in zip(stats[ii],precision):
		sl=sl+'&'		
		sl=sl+sstr(st,jj)
        sl=sl+'\\\\'
        s.append(sl)
s.append( ' \end{tabular} ' )
s.append( ' \end{table} ' )


f = open("./tables_tex/table"+str(vec[0])+'_'+str(vec[-1])+'.txt', "w")
for st in s:
	print st
        f.write(st+'\n')
f.close()









from netCDF4 import Dataset as NF
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as matplotlib
import numpy as np
import scipy.interpolate as sci_interp
import mod_data as md
#import gshhs as gshhs
from matplotlib.lines import Line2D

ddir='./';
grdname=ddir+'grd2.nc';

lonmin =-11;
lonmax = 1;
latmin = 34.7;
latmax = 38.5;
lonmin =-7.05;
lonmax = 2.5;
lonmax = 5.2;
latmin = 35.05;
latmax = 36.9;

# Isotropic grid
gibRadEarth=6371e3;
#gibDx=4e3;
gibDx=2e3;
dl = ((gibDx/gibRadEarth)/np.cos(36*np.pi/180))*180/np.pi;

lon=np.arange(lonmin,lonmax,dl);

i=0;
lat=[latmin];
while lat[-1]<=latmax:
    i=i+1;
    b=[lat[i-1]+dl*np.cos(lat[i-1]*np.pi/180)]
    lat=np.r_[lat,b];

ff = NF('./ETOPO1_Ice_g_gdal.grd', 'r')

x=np.linspace(-180.,180.,360*60+1);
y=np.linspace(-90.,90.,180*60+1);

xi=np.where((x>lonmin) & (x < lonmax))[0];
xi=np.arange(xi[0]-4,xi[-1]+4)
yi=np.where((y>latmin) & (y < latmax))[0];
yi=np.arange(yi[0]-4,yi[-1]+4)

x1=np.arange(np.size(yi))
x2=np.arange(np.size(yi))

for j in np.arange(0,np.size(yi)):
    x1[j]=np.size(x)*(np.size(y)-yi[j])+xi[0];
    x2[j]=np.size(x)*(np.size(y)-yi[j])+xi[-1];

tz=np.nan*np.ones([np.size(yi),np.size(xi)]);

if 1:
    for j in np.arange(0,np.size(yi)):
        tz[j]=ff.variables['z'][x1[j]:x2[j]+1];

f = sci_interp.RectBivariateSpline(x[xi],y[yi], tz.transpose())
h = -f(lon, lat).transpose()

md.create_grd(grdname,lat,lon,'l')

grid=NF(grdname,'a')
grid.variables['lonh'][:]=lon
grid.variables['lath'][:]=lat
grid.variables['D'][:]=h
grid.close()

print "lenlat: " + str(lat[-1]-lat[0])
print "lenlon: " + str(lon[-1]-lon[0])
print ""
print "lowlat: " + str(min(lat))
print "westlon: " + str(min(lon))

print "size: " + str(h.shape)

dy=(np.diff(lat)*np.pi/180)*gibRadEarth

print "dy min: "+str( np.min(dy) )
print "dy max: "+str( np.max(dy) )

axlo=x[xi];
axla=y[yi];

plt.close('all')

low1=0.6
low2=0.2
h1=0.375

fac=6
fig=plt.figure(figsize=(fac*2+2,fac*1))
of1=0.32-0.32*12./(14)
ax1 = plt.axes([0.045,low1,0.54+of1,h1])
plt.figtext(0.02,low1+h1,'a)',fontsize=15)
#ax1=plt.subplot(211)
h=-h
h[h>=0.]=10.0
norm1=matplotlib.colors.Normalize(vmin = -1500, vmax = 0.1, clip = False)
palette = matplotlib.cm.jet
palette.set_over('gray', 0.5)

cf=plt.contourf(lon,lat,h,\
    np.linspace(-1500,0.,31),\
    cmap=palette,\
    norm = norm1, extend='both'
    )
ct=plt.contour(lon,lat,h,[-100.,-500.,-1000.],colors='k',linestyles='solid')
plt.xlim(lon[0], lon[-1])
plt.ylim(lat[0], lat[-1])
#cb=plt.colorbar(cf)
#cb.set_label('Depth [m]')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

los=17
loe=los+517
las=9
lae=las+85
plt.hlines([lat[las],lat[lae]],lon[los],lon[loe],linestyles='--',linewidth=2)
plt.vlines([lon[los],lon[loe]],lat[las],lat[lae],linestyles='--',linewidth=2)

tr1=ax1.transData.transform((lon[los],lat[las]))
tr2=ax1.transData.transform((lon[loe],lat[lae]))
xx1=tr1[0]; xx2=tr2[0]
inv = fig.transFigure.inverted()
tr1=inv.transform(tr1)
tr2=inv.transform(tr2)

############################################
gname='../run170/data/grd.nc';
ff=NF(gname,'r')
xh=ff.variables['xh'][:]
yh=ff.variables['yh'][:]
D=ff.variables['D'][:]
ff.close()
Ds=284.
Dn=880.

isill=np.argwhere(D[np.argmax(D[:,0]),:]==Ds)[0]
inarrow=np.argwhere(D[np.argmax(D[:,0]),:]==Dn)[0]

D=-D
D[D>=-0.1]=10.

#ax2=plt.subplot(212)
ax2 = plt.axes([0.07,low2,0.5+of1,0.4])

pts=plt.gca().get_position().get_points()
bb_new=matplotlib.transforms.Bbox([  [tr1[0],pts[0,1]]  ,  [tr2[0],pts[0,1]+tr2[1]-tr1[1]] ])
plt.gca().set_position(bb_new)
cf=plt.contourf(xh*1e-3,yh*1e-3,D,\
    np.linspace(-1500,0.,31),\
    cmap=palette,\
    norm = norm1, extend='both'
    )
ct=plt.contour(xh*1e-3,yh*1e-3,D,[-100.,-500.,-1000.],colors='k',linestyles='solid')

lonsponge=[plt.xlim()[0]+30.,  plt.xlim()[1]-300.]
invax=ax2.transAxes.inverted()
trsp1=ax2.transData.transform( (lonsponge[1],0.) )
trsp2=invax.transform(trsp1)
plt.vlines(np.linspace(trsp2[0],1,20),0,1,transform=ax2.transAxes,\
       linestyles='dotted',linewidth=2,color='w')
trsp1=ax2.transData.transform( (lonsponge[0],0.) )
trsp2=invax.transform(trsp1)
plt.vlines(np.linspace(0,trsp2[0],4),0,1,transform=ax2.transAxes,\
       linestyles='dotted',linewidth=2,color='w')

plt.xlabel('x [km]')
plt.ylabel('y [km]') 

ax3 = plt.axes([0,0,1,1], axisbg=(1,1,1,0),frame_on=False)
tr4=ax1.transData.transform( (-4,36.3) )
tr5=ax2.transData.transform( (xh[isill]*1e-3,0) )
tr3=ax2.transData.transform( (xh[inarrow]*1e-3,0) )
tr4=inv.transform(tr4)
tr5=inv.transform(tr5)
tr3=inv.transform(tr3)
#plt.vlines(tr5[0],tr5[1],tr4[1],transform=ax3.transAxes)
#plt.vlines(tr3[0],tr5[1],tr4[1],transform=ax3.transAxes)



ax3.xaxis.set_ticks_position("none")
ax3.yaxis.set_ticks_position("none")

p1,p2,p3,p4=ax2.get_position().bounds
cax=fig.add_axes( [p1,0.08,p3,0.1*p4] )
myformat = plt.FormatStrFormatter('%2.0f')
cb=plt.colorbar(cf,cax,orientation='horizontal',\
    format=myformat,ticks=np.linspace(-1500,0.,7))
cb.set_label('Depth [m]')

##############################################
ax4 = plt.axes([0.69,low1,0.35-of1,h1])
plt.figtext(0.645,low1+h1,'c)',fontsize=15)
norm1=matplotlib.colors.Normalize(vmin = -1000, vmax = 0.1, clip = False)
palette = matplotlib.cm.jet
palette.set_over('gray', 0.5)

sx=np.arange(42,92)
sx=np.arange(32,102)
sy=np.arange(37,68)
if 1:
    cf=plt.contourf(lon[sx],lat[sy],h[sy,sx[0]:sx[-1]+1],\
        np.linspace(-1000,0.,21),\
        cmap=palette,\
        norm = norm1, extend='both'
        )
ct=plt.contour(lon[sx],lat[sy],h[sy,sx[0]:sx[-1]+1],[-100.,-500.,-1000.],colors='k',linestyles='solid')

#window=( [lon[sx[0]],lon[sx[-1]],lat[sy[0]],lat[sy[-1]]] )
#mp = gshhs.gshhs(window=window, resolution='h', clip=False, area_thresh=0., max_level=1)
   
#for lon, lat, level in zip(mp.lon, mp.lat, mp.level):        
#    ax4.add_line(Line2D(lon, lat,color='k')) 

plt.ylabel('Latitude')
plt.grid(True)

################################################
ax5 = plt.axes([0.69,low2,0.35-of1,tr2[1]-tr1[1]])
plt.figtext(0.02,low2+tr2[1]-tr1[1],'b)',fontsize=15)
plt.figtext(0.64,low2+tr2[1]-tr1[1],'d)',fontsize=15)
lx=np.arange(16,16+len(sx))
ly=np.arange(27,27+len(sy))
if 1:
    cf=plt.contourf(xh[lx]*1e-3,yh[ly]*1e-3,D[ly,lx[0]:lx[-1]+1],\
        np.linspace(-1000,0.,21),\
        cmap=palette,\
        norm = norm1, extend='both'
        )
    ct=plt.contour(xh[lx]*1e-3,yh[ly]*1e-3,D[ly,lx[0]:lx[-1]+1],[-100.,-500.,-1000.],colors='k',linestyles='solid')

else:
    extent = (xh[lx[0]]*1e-3-1,xh[lx[-1]]*1e-3+1,yh[ly[0]]*1e-3-1,yh[ly[-1]]*1e-3+1)
    cf=plt.imshow(D[ly,lx[0]:lx[-1]+1],interpolation='nearest',extent=extent,\
    cmap=palette,\
    norm = norm1)
plt.xlabel('x [km]')
plt.ylabel('y [km]') 
plt.yticks(np.linspace(-20,20,5))
plt.grid(True)
plt.annotate('Sill', [xh[isill]*1e-3,0], xytext=[xh[isill]*1e-3,-18], xycoords='data',
         arrowprops=dict(arrowstyle="->",facecolor='black'),textcoords='data',horizontalalignment='center')
plt.annotate('Contraction', [xh[inarrow]*1e-3,-8], xytext=[xh[inarrow]*1e-3,-28], xycoords='data',
         arrowprops=dict(arrowstyle="->",facecolor='black'),textcoords='data',horizontalalignment='center')
plt.annotate('Exit', [xh[inarrow+9]*1e-3,12], xytext=[xh[inarrow+4]*1e-3,21], xycoords='data',
         arrowprops=dict(arrowstyle="->",facecolor='black'),textcoords='data',horizontalalignment='center')

ax6 = fig.add_axes(ax3.get_position(), frameon=False)
tr6=ax4.transAxes.transform( (0.,1.) )
tr5=ax5.transData.transform( (xh[isill]*1e-3,0) )
tr3=ax5.transData.transform( (xh[inarrow]*1e-3,0) )
tr99=ax5.transData.transform( (xh[inarrow+9]*1e-3,0) )
tr6=inv.transform(tr6)
tr5=inv.transform(tr5)
tr3=inv.transform(tr3)
tr99=inv.transform(tr99)
plt.vlines(tr5[0],tr5[1],tr6[1],transform=ax6.transAxes,linestyles='--',linewidth=2)
plt.vlines(tr3[0],tr5[1],tr6[1],transform=ax6.transAxes,linestyles='--',linewidth=2)
plt.vlines(tr99[0],tr5[1],tr6[1],transform=ax6.transAxes,linestyles='--',linewidth=2)

ax6.xaxis.set_ticks_position("none")
ax6.yaxis.set_ticks_position("none")
p1,p2,p3,p4=ax5.get_position().bounds
cax=fig.add_axes( [p1,0.08,p3,0.1*p4] )
myformat = plt.FormatStrFormatter('%2.0f')
cb=plt.colorbar(cf,cax,orientation='horizontal',format=myformat,\
    ticks=np.linspace(-1000,0.,6),extend='both')
cb.set_label('Depth [m]')

xt=ax4.axes.get_xticks()
pad=plt.rcParams['xtick.major.pad']
xtl=ax4.axes.set_xticklabels([])

for ii in range(1,len(xt)-1):
     px=ax4.transData.transform((xt[ii],0))
     py=ax4.transAxes.transform((xt[ii],0))
     p=[px[0],py[1]-pad] #display coordinates 
     p2=inv.transform(p) #figure coordinates 
     plt.figtext(p2[0],p2[1],xt[ii],ha='center',va='top',backgroundcolor='white')
py=ax4.transAxes.transform((0.5,0))
p3=inv.transform((py[0],0))
plt.figtext(p3[0],p2[1]-0.04,'Longitude',ha='center',va='top',backgroundcolor='white')

lat36=np.argmin(np.abs(lat-36.))

plt.savefig('./figures/fig1.pdf')

from netCDF4 import Dataset as NF
import matplotlib.pyplot as plt
import numpy as np
import mod_data as md
import mod_triang as tr
import os as os
import inspect as inspect

rundir= os.path.dirname(inspect.getfile(inspect.currentframe()))
grdname=rundir+'/../data/grd.nc';


#set up basin
domx=750e3 #distance Cam. sill-Gibraltar ~40km
domx+=-10.723*(15+2*180)
domy=170e3
#domy+=2e3
depth=2e3

north=tr.make_coast_rect()
bottom=tr.make_bottom()
print north.shape
print bottom.shape
#res_c=15
#res_b=180.
#nx=res_c+2.*(res_b-2.);
nx=float(north.shape[0])
ny=np.round((domy/domx)*nx)
#
xh=domx*(np.arange(nx)/nx-0.5)
yh=domy*(np.arange(ny)/ny)
yh=yh-yh[np.floor(yh.shape[0]/2.)]
#

narrow_i=np.argmin(north)

ex=np.argwhere(np.diff(north)>60e3)[0][0] #strait exit

bot=np.nan*np.ones((ny,nx))
for jj in np.arange(bot.shape[0]):
    for ii in np.arange(ex+1):
        dist=north[ii]-np.abs(yh[jj])
#+0.5*np.diff(yh)[0]
        if dist<0.:
            bot[jj,ii]=0.
        else:
            bot[jj,ii]=bottom[ii]*(1.-np.abs(((dist-north[ii])/north[ii])))

        
mid_i=np.argmax(bot[:,0]).flatten()
#gr=np.sort(bot[:,ex]); d1=gr[gr!=0.][0]
d1=20.

d2=100.
d3=500.
d4=920.
dxd2=4 #distance coast to d2 isobath in grdpts
height1=d1+np.linspace(0.,1.,dxd2)*(d2-d1)
dxd3=14 #d3 isobath 
height2=d2+np.linspace(0.,1.,dxd3)*(d3-d2)
dxd4=4 
height3=d3+np.linspace(0.,1.,dxd4)*(d4-d3)
height=np.r_[height1[:-1],height2[:-1],height3]
w=len(height) #width
#height=50.+height
i1=mid_i-5 - np.round(np.linspace(0.,1.,w)*(15))
i2=mid_i
for ii in range(w):
    bot[:i1[ii],ex+1+ii]=height[ii]
    bot[bot.shape[0]-i1[ii]:,ex+1+ii]=height[ii]

    zo=np.linspace(0.,1.,i2-i1[ii]+1)
    #if ii>0:
    zo=zo**3
    tmp=height[ii]+zo*(height[-1]-height[ii])
    bot[i1[ii]:i2+1,ex+1+ii]=tmp
    bot[i2:bot.shape[0]-i1[ii],ex+1+ii]=tmp[-1::-1]


for jj in np.arange(bot.shape[0]):
    for ii in np.arange(ex+w,bot.shape[1]):
        bot[jj,ii]=bottom[ii-w]

# fix zo
for ii in range(w):
    dif=bot[:,ex+1+ii]-bot[:,ex+ii]
    bot[:,ex+1+ii][dif<0]=bot[:,ex+ii][dif<0]
       
sill_i=np.argmin(bottom)

plt.close()
plt.close()
plt.close()

plt.close()
plt.close()
plt.close()
plt.figure()
plt.plot(north)
plt.axvline(sill_i, ymin=0, ymax=1,color='k')
plt.axvline(narrow_i, ymin=0, ymax=1,color='k')
plt.axvline(ex, ymin=0, ymax=1,color='r')
plt.figure()
plt.plot(bottom,'.')
plt.axvline(sill_i, ymin=0, ymax=1,color='k')
plt.axvline(narrow_i, ymin=0, ymax=1,color='k')

plt.figure()
#bot[bot>500]=np.nan
plt.imshow(bot,interpolation='nearest',origin='bottom')
plt.colorbar()

print "Depth at Sill:    " +str(bottom[sill_i])
print "Depth at Narrows: " +str(bottom[narrow_i])
W=2*north
print "Width at Sill:    " +str(W[sill_i])
print "Width at Narrows: " +str(W[narrow_i])

print "Distance Sill-Narrow: " +str((domx/nx)*(narrow_i-sill_i))

print "dx (domx/nx): "+str(domx/nx)
print "dy (domy/ny): "+str(domy/ny)


# run90+
dxh=np.diff(xh)
xh=np.r_[xh,xh[-1]+np.cumsum(dxh[-100:])]
bot=np.c_[bot,bot[:,-100:]]
xh=np.r_[xh,xh[-1]+np.cumsum(dxh[-60:])]
bot=np.c_[bot,bot[:,-60:]]
# run123+

xh=np.r_[xh,xh[-1]+np.cumsum(dxh[-100:])]
bot=np.c_[bot,bot[:,-100:]]
#run135+
xh=np.r_[xh,xh[-1]+np.cumsum(dxh[-24:])]
bot=np.c_[bot,bot[:,-24:]]
    
if 0: #large domain
    xs=140.
    xh=xh[xs:]
    bot=bot[:,xs:]
if 0: #small
    xs=150.; xe=xs+90.
    xh=xh[xs:xe]
    bot=bot[:,xs:xe]
if 0: #channel
    xs=160.; xe=xs+70.
    w=10.
    ys=43.-w; ye=43.+w-1;
    xh=xh[xs:xe]
    yh=yh[ys:ye]
    bot=bot[ys:ye,xs:xe] 

if 0: #Alboran
    xs=160.; 
    xh=xh[xs:]
    bot=bot[:,xs:]   
    
if 1: #Alboran run90+
    xs=140.; 
    xh=xh[xs:]
    bot=bot[:,xs:]

bot[bot==0.]=0.1
bot[0,:]=0.1
bot[-1,:]=0.1 # dirty 

###############################
# shelf on southern bdy
Dn=880.
jmid=np.argmax(bot[:,0])
inarrow=np.argwhere(bot[jmid,:]==Dn)[0][0]
iexit=np.argwhere(bot[10,inarrow:]>1)[0][0]+inarrow-1 #index of strait exit
i1500=np.argwhere(bot[jmid,:]==1500)[0][0]
s1=bot[2,iexit:i1500].copy() 

nn=len(s1)
##[sq,tr]=np.meshgrid(s1,np.arange(nn))

s2=[[s1[ii]]* (bot.shape[1]-i1500) for ii in range(len(s1))]
s2=np.array(s2)
bot[:s2.shape[0],i1500:bot.shape[1]]=s2

sq=bot[:nn,iexit:i1500]
co=0
for i in np.arange(nn):
    for j in np.arange(1+co,nn):
        sq[i,j]=s1[co]
    co=co+1
#bot[:s1.shape[0],iexit:iexit+nn]=sq
###############################
Ds=284.
isill=np.argwhere(bot[jmid,:]==Ds)[0]
#tres forcas
itf=isill+124
lj=[13,13,12,11,9,7,5,5,4,4,3,3,2,2]
for ii in range(len(lj)):
    gr=50
    tmp=bot[:gr,itf+ii]
    bot[:lj[ii],itf+ii]=0.1
    bot[:lj[ii],itf-ii]=0.1
    #bot[lj[ii]-1:lj[ii]+gr-1,itf+ii]=tmp #unstable
    #bot[lj[ii]-1:lj[ii]+gr-1,itf-ii]=tmp
    bot[lj[ii]:lj[ii]+gr-2,itf+ii]=tmp[2:]
    bot[lj[ii]:lj[ii]+gr-2,itf-ii]=tmp[2:]    



## cape_v2 smooth shelf around cape
#sl=np.linspace(0,1,20)
#for ii in range(1,13):
#    #w2=np.where(bot[ii,-1::-1]==0.1)[0][0]
#    ib=itf-len(lj)-2
#    w=np.where(bot[ii,ib:]==0.1)[0]+ib
#    w0=w[0]
#    w1=w[-1]
#    gu=np.r_[sl,np.ones(w1-w0),sl[-1::-1]]
#    de=bot[ii,w0-len(sl)]
#    bot[ii,w0-len(sl):w1+len(sl)]=0.1*gu+de*(1.-gu)

#for ii in range(13,30):
#    de=bot[ii,w0-len(sl)]
#    nu=(ii-12.)/(30-13)
#    print nu
#    de2=de*nu+0.1*(1.-nu)
#    de2=bot[29,w0]*nu+0.1*(1.-nu)
#    bot[ii,w0-len(sl):w1+len(sl)]=de2*gu+de*(1.-gu)

#de2=bot[0:30,-1]
#gr=de2.shape[0]
#for ii in range(gr):
#    de=bot[ii+13,w0-len(sl)]
#    bot[ii+13,w0-len(sl):w1+len(sl)]=de2[ii]*gu+de*(1.-gu)




md.create_grd(grdname,yh,xh,'x')

grid=NF(grdname,'a')
grid.variables['xh'][:]=xh
grid.variables['yh'][:]=yh
grid.variables['D'][:]=bot
grid.close()
#
print ""
print "lenlat: " + str(np.diff(yh)[0]*yh.shape[0])
print "lenlon: " + str(np.diff(xh)[0]*xh.shape[0])
print ""
print "lowlat: " + str(min(yh))
print "westlon: " + str(min(xh))
#
print "size: " + str(bot.shape)



plt.figure()
bot[bot==0.1]=np.nan
#plt.imshow(bot,interpolation='nearest')
plt.contourf(xh*1e-3,yh*1e-3,bot,30)
plt.colorbar()

plt.figure()
plt.imshow(bot[-1::-1,:])





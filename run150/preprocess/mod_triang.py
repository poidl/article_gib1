import numpy as np

#def make_bottom():
#    Ds=284.
#    Dn=880.
#    depth=2e3
#    res_c=15.
#    res_b=180.    
#    res=res_c+2.*res_b    
#    
#    sig_sill=0.5 
#    narrow_i=int(0.99*res_c); narrow_i+=res_b-1
#    sill_i=int(0.3*res_c); sill_i+=res_b-1
#    
#    xdist=np.linspace(0.,1.,res)
#    gau=np.exp(-(xdist-xdist[sill_i])**2/(2*sig_sill**2))
#    c=(Ds-Dn)/(gau[sill_i]-gau[narrow_i])
#    D=Dn-c*gau[narrow_i]+c*gau   
#    D[D>depth]=depth
#    
#    # Mediterranean
#    ind=res_b-3
#    vec=np.linspace(0.,1.,ind)
#    tmp=D[-ind]+(depth-D[-ind])*vec
#    D=np.r_[D[:-ind], tmp]
#    #    #smooth
#    #for ii in np.arange(-2.,3.):    
#    #    D[-ind+ii]=0.5*(D[-ind+ii]+D[-ind+ii-1])
#    
#    # Atlantic
#    ind=res_b
#    vec=np.linspace(0.,1.,ind)
#    tmp=D[ind-1]+(depth-D[ind-1])*vec[-1::-1]
#    D=np.r_[tmp,D[ind:]]
#    
#    
#    return D
    


def make_bottom():
    # like make_bottom, but make slope in Atlantic steeper, i.e.
    Ds=284.
    Dn=880.
    depth=2e3
    res_c=15.
    res_b=180.    
    res=res_c+2.*res_b    
    
    sig_sill=0.5 
    narrow_i=int(0.99*res_c); narrow_i+=res_b-1
    sill_i=int(0.3*res_c); sill_i+=res_b-1
    
    xdist=np.linspace(0.,1.,res)
    gau=np.exp(-(xdist-xdist[sill_i])**2/(2*sig_sill**2))
    c=(Ds-Dn)/(gau[sill_i]-gau[narrow_i])
    D=Dn-c*gau[narrow_i]+c*gau   
    D[D>depth]=depth
    
    # Mediterranean
    ind=res_b+1
    vec=np.linspace(0.,1.,ind)
    tmp=920.+(6e3-920.)*vec
    tmp[tmp>1500.]=1500.
    flat=9
    t2=np.repeat(920.,flat-1)
    t2=np.r_[D[-ind-1],t2]
    D=np.r_[D[:-ind-1], t2,tmp[:-flat+1]]
    
    #    #smooth
    #for ii in np.arange(-2.,3.):    
    #    D[-ind+ii]=0.5*(D[-ind+ii]+D[-ind+ii-1])
    
    # Atlantic
    ind=res_b
    vec=np.linspace(0.,1.,ind)
    tmp=D[ind-1]+(0.5*depth-D[ind-1])*vec[-1::-1]
    D=np.r_[tmp,D[ind:]]
    
    
    return D    
    
    
#    depth=2e3
#    xdist=np.linspace(0.,1.,res_b)
#    gau1=xdist
#    db=-(D[-1]-D[-2])/np.diff(xdist)[0]
#   
#    b=2.*(D[-1]-depth-0.5*db)
#    a=0.5*(db-b)
#    c=depth
#    gau2=a*xdist[-1::-1]**2 + b*xdist[-1::-1] + c
#
#    gau2=D[-1]+(depth-D[-1])*gau1
#
#
#    bottom=np.r_[D[0]+(depth-D[0])*gau1[-1::-1],D[1:-1], gau2]
#    


def make_coast2():
    Ws=22.3e3
    Wn=13.8e3   
    
    domy=170e3
    
    res_c=15.

    sig_narrow=0.2
    narrow_i=int(0.99*res_c)
    sill_i=int(0.3*res_c)
    
    xdist=np.linspace(0.,1.,res_c)
    gau=np.exp(-(xdist-xdist[narrow_i])**2/(2*sig_narrow**2))
    c=(Ws-Wn)/(gau[sill_i]-gau[narrow_i])
    W=Ws+c*gau-c*gau[sill_i]

    
    res_b=180.
    sig_basin=0.08
    xdist=np.linspace(0.,1.,res_b)
    gau=1.-np.exp(-(xdist)**2/(2*sig_basin**2))
    north=0.5*np.r_[W[0]+(domy-W[0])*gau[-1::-1],W[1:-1], W[-1]+(domy-W[-1])*gau]

    return north
    
def make_coast_rect():
    Ws=22.3e3
    Wn=13.8e3   
    
    domy=170e3
    
    res_c=15.

    sig_narrow=0.2
    narrow_i=int(0.99*res_c)
    sill_i=int(0.3*res_c)
    
    xdist=np.linspace(0.,1.,res_c)
    gau=np.exp(-(xdist-xdist[narrow_i])**2/(2*sig_narrow**2))
    c=(Ws-Wn)/(gau[sill_i]-gau[narrow_i])
    W=Ws+c*gau-c*gau[sill_i]

    ca=4
    W=np.r_[W,W[-ca:][-1::-1]]
    
    
    res_b=180.
    sig_basin=0.08
    xdist=np.linspace(0.,1.,res_b)
    gau=1.-np.exp(-(xdist)**2/(2*sig_basin**2))
    cl=6 #contraction length (in gridpts)
    north=0.5*np.r_[W[0]+(domy-W[0])*gau[-1::-1],W[1:-1],np.repeat(W[-1],cl), W[-1]+(domy-W[-1])*np.ones(res_b-cl-ca)]

    return north
    
    
    
def make_bottom2():
    Ds=284.
    Dn=880.    
    res_c=15.
    
    sig_sill=0.8
    narrow_i=int(0.99*res_c)
    sill_i=int(0.3*res_c)
    
    xdist=np.linspace(0.,1.,res_c)
    gau=np.exp(-(xdist-xdist[sill_i])**2/(2*sig_sill**2))
    c=(Ds-Dn)/(gau[sill_i]-gau[narrow_i])
    D=Dn-c*gau[narrow_i]+c*gau   
    
    
    res_b=180.
    depth=2e3
    xdist=np.linspace(0.,1.,res_b)
    gau1=xdist
    db=-(D[-1]-D[-2])/np.diff(xdist)[0]
   
    b=2.*(D[-1]-depth-0.5*db)
    a=0.5*(db-b)
    c=depth
    gau2=a*xdist[-1::-1]**2 + b*xdist[-1::-1] + c
    
    gr=920.
    gau2=D[-1]+(depth-D[-1])*gau1


    bottom=np.r_[D[0]+(depth-D[0])*gau1[-1::-1],D[1:-1], gau2]

#    Ds=284.
#    Dn=880.    
#    res_c=15.
#    sig_sill=0.8
#    narrow_i=int(0.99*res_c)
#    sill_i=int(0.3*res_c)
#    
#    xdist=np.linspace(0.,1.,res_c)
#    gau=np.exp(-(xdist-xdist[sill_i])**2/(2*sig_sill**2))
#    c=(Ds-Dn)/(gau[sill_i]-gau[narrow_i])
#    D=Dn-c*gau[narrow_i]+c*gau   
#    
#    
#    res_b=180.
#    depth=2e3
#    xdist=np.linspace(0.,1.,res_b)
#    gau1=xdist
#    db=-(D[-1]-D[-2])/np.diff(xdist)[0]
#   
#    b=2.*(D[-1]-depth-0.5*db)
#    a=0.5*(db-b)
#    c=depth
#    gau2=a*xdist[-1::-1]**2 + b*xdist[-1::-1] + c
#
#    gau2=D[-1]+(depth-D[-1])*gau1
#
#
#    bottom=np.r_[D[0]+(depth-D[0])*gau1[-1::-1],D[1:-1], gau2]

    return bottom
    
    

def make_channel(res_c,sig_narrow,narrow_i,sig_sill,sill_i):
    Ds=284.
    Dn=880.
    Ws=22.3e3
    Wn=13.8e3
    
    xdist=np.linspace(0.,1.,res_c)
     
    #width: gau is gaussian and exactly 0.0 at sill_i and 1.0 at narrow_i
    gau=np.exp(-(xdist-xdist[narrow_i])**2/(2*sig_narrow**2))
    c=(Ws-Wn)/(gau[sill_i]-gau[narrow_i])
    W=Ws+c*gau
    
    #bottom
    gau=np.exp(-(xdist-xdist[sill_i])**2/(2*sig_sill**2))
    c=(Ds-Dn)/(gau[sill_i]-gau[narrow_i])
    D=Dn-c*gau[narrow_i]+c*gau
    
    return (W,D)


def make_basin(res_b,depth,W,D,sig_basin,domy):
    
    xdist=np.linspace(0.,1.,res_b)
    #gau3 is exactly zero at eastern bdy
    gau=1.-np.exp(-(xdist)**2/(2*sig_basin**2))
    north=0.5*np.r_[W[0]+(domy-W[0])*gau[-1::-1],W[1:-1], W[-1]+(domy-W[-1])*gau]

    
    gau1=xdist
    #gau2=1.-np.exp(-(xdist)**2/(2*sig_basin**2))
    db=-(D[-1]-D[-2])/np.diff(xdist)[0]
    #bottom=np.r_[D[0]+(depth-D[0])*gau1[-1::-1],D[1:-1], D[-1]+(depth-D[-1])*gau2]
    
    b=2.*(D[-1]-depth-0.5*db)
    a=0.5*(db-b)
    c=depth
    gau2=a*xdist[-1::-1]**2 + b*xdist[-1::-1] + c

    gau2=D[-1]+(depth-D[-1])*gau1

    #ii=20
    #ls=np.linspace(0.,1.,res_b-ii)
    #gau2[ii:]=gau2[ii]+(depth-gau2[ii])*ls
    bottom=np.r_[D[0]+(depth-D[0])*gau1[-1::-1],D[1:-1], gau2]
    return (north,bottom)
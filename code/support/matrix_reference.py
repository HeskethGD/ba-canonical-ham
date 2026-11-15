"""Independent numerical matrix reference; rates use geometric time."""
import numpy as np
s2=np.sqrt(2)

def build(L,M,N,bL,bM,bN,u,us,v,vs,w,ws):
    # dagger = transpose of starred array
    ud,vd,wd = us.T, vs.T, ws.T
    X=np.block([[-2*np.diag(bL), s2*u, wd],[s2*ud, -2*np.diag(bM), v],[-w, -vd, np.diag(bN)]])
    D=np.diag(np.r_[np.ones(L),-np.ones(M),np.zeros(N)])
    return X,D

def T(L,M,N,bL,bM,v,vs,w,ws,nu):
    return np.block([[-np.diag(bL), np.zeros((L,M)), ws.T],[np.zeros((M,L)), np.diag(bM), -v],[-w, vs.T, -nu*np.eye(N)]])

def rates(L,M,N,bL,bM,bN,u,us,v,vs,w,ws):
    BL,BM,BN=np.diag(bL),np.diag(bM),np.diag(bN)
    ud,vd,wd = us.T, vs.T, ws.T
    du = BL@u+u@BM+s2*wd@vd
    dud= -BM@ud-ud@BL-s2*v@w
    dv = BM@v+v@BN+s2*ud@wd
    dvd= -BN@vd-vd@BM-s2*w@u
    dw = BN@w+w@BL+s2*vd@ud
    dwd= -BL@wd-wd@BN-s2*u@v
    return du,dud.T,dv,dvd.T,dw,dwd.T

def H(L,M,N,bL,bM,bN,u,us,v,vs,w,ws):
    h=0
    for l in range(L):
        for m in range(M): h+=(bL[l]+bM[m])*u[l,m]*us[l,m]
    for m in range(M):
        for n in range(N): h+=(bM[m]+bN[n])*v[m,n]*vs[m,n]
    for n in range(N):
        for l in range(L): h+=(bN[n]+bL[l])*w[n,l]*ws[n,l]
    h+=s2*(np.trace(u@v@w)+np.einsum('lm,mn,nl',us,vs,ws))
    return h

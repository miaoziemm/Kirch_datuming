#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validate phase-residual butterfly for three current wavefield extrapolation kernels:
1. receiver-side full kernel
2. source-side cmp=0 full kernel
3. source-side cmp=1 common-receiver regrouped source-axis kernel

This script is self-contained. It mimics the MATLAB phase-residual butterfly structure,
but replaces the MATLAB test kernel ctheta*exp(+i*w*t) with the current extrapolation
kernel K = W * AA * H_fast.
"""
import math
import numpy as np
from dataclasses import dataclass

def relerr(a,b):
    return np.linalg.norm(a-b)/(np.linalg.norm(b)+1e-30)

def cheb_nodes(order,a,b):
    c=0.5*(a+b); h=0.5*(b-a)
    return np.array([c+h*math.cos(math.pi*(2*k+1)/(2*order)) for k in range(order)],float)

def bary_weights(nodes):
    nodes=np.asarray(nodes,float)
    q=len(nodes)
    w=np.ones(q,float)
    for j in range(q):
        prod=1.0
        for k in range(q):
            if k!=j:
                prod *= nodes[j]-nodes[k]
        w[j]=1.0/prod
    return w

def lagrange_all(nodes, weights, x):
    nodes=np.asarray(nodes,float); weights=np.asarray(weights,float)
    diff=x-nodes
    hit=np.where(np.abs(diff)<1e-9)[0]
    if len(hit):
        L=np.zeros(len(nodes)); L[hit[0]]=1.0; return L
    tmp=weights/diff
    return tmp/tmp.sum()

@dataclass
class Box:
    id:int; level:int; begin:int; end:int; parent:int; left:int; right:int
    xmin:float; xmax:float; center_index:int; center_coord:float
    cheb_index:np.ndarray; cheb_coord:np.ndarray; bary_w:np.ndarray

class Tree1D:
    def __init__(self,n,x0,dx,leaf_size,order):
        self.n=n; self.x0=x0; self.dx=dx; self.leaf_size=leaf_size; self.order=order
        self.boxes=[]; self.depth=0
        self._build(0,n,0,-1)

    def _assign_cheb(self,begin,end):
        xmin=self.x0+begin*self.dx
        xmax=self.x0+(end-1)*self.dx
        nbox=end-begin
        q=min(self.order,nbox)
        raw=cheb_nodes(q,xmin,xmax)
        used=set(); inds=[]; coords=[]
        for rn in raw:
            idx=int(math.floor((rn-self.x0)/self.dx+0.5))
            idx=max(begin,min(end-1,idx))
            if idx in used:
                best=None; bestd=10**9
                for j in range(begin,end):
                    if j not in used:
                        d=abs(j-idx)
                        if d<bestd:
                            best=j; bestd=d
                idx=best if best is not None else idx
            used.add(idx)
            inds.append(idx); coords.append(self.x0+idx*self.dx)
        coords=np.asarray(coords,float)
        return np.asarray(inds,int), coords, bary_weights(coords)

    def _build(self,begin,end,level,parent):
        bid=len(self.boxes)
        xmin=self.x0+begin*self.dx
        xmax=self.x0+(end-1)*self.dx
        center_index=int(round((begin+end-1)/2))
        center_coord=self.x0+center_index*self.dx
        ci,cc,bw=self._assign_cheb(begin,end)
        box=Box(bid,level,begin,end,parent,-1,-1,xmin,xmax,center_index,center_coord,ci,cc,bw)
        self.boxes.append(box)
        self.depth=max(self.depth,level)
        if end-begin>self.leaf_size:
            mid=begin+(end-begin)//2
            l=self._build(begin,mid,level+1,bid)
            r=self._build(mid,end,level+1,bid)
            self.boxes[bid].left=l; self.boxes[bid].right=r
        return bid

    def level_boxes(self, level):
        return [b.id for b in self.boxes if b.level==level]
    def root(self):
        return 0

def build_filter_F(tau,dt,nsam):
    G=np.zeros(nsam,float)
    for k in range(1,nsam):
        x=(tau+k*dt)/tau
        G[k]=math.sqrt(x*x-1.0)
    for k in range(nsam-1):
        G[k]=G[k+1]-G[k]
    for k in range(nsam-2,0,-1):
        G[k]=G[k]-G[k-1]
    return G[:nsam-1].copy()

def Hfast(tau,iw,nfft,dt,nsam):
    if tau<=0 or not np.isfinite(tau):
        return 0j
    F=build_filter_F(tau,dt,nsam)
    q=tau/dt
    m=math.ceil(q)
    delta=m-q
    omega=2*math.pi*iw/nfft
    k=np.arange(nsam-1)
    sumF=np.sum((F/dt)*np.exp(-1j*omega*k))
    interp=(1-delta)*np.exp(-1j*omega*m)+delta*np.exp(-1j*omega*(m-1))
    return interp*sumF

def phase_ratio(Ha,Hb):
    z=Ha*np.conj(Hb)
    mag=abs(z)
    if mag<1e-20:
        return 1+0j
    return z/mag

def taper_weight(left,center,right,tap,stride=1):
    if tap<=0:
        return 1.0
    # Preserve the integer-division semantics discussed for the C/C++ code.
    wl=1.0 if (center-left)>=tap else float(((center-left)//stride)//tap)
    wr=1.0 if (right-center)>=tap else float(((right-center)//stride)//tap)
    return wl*wr

def table_slope_first_index(table,row,col,dcoord):
    n=table.shape[0]
    if n<=1 or dcoord==0:
        return 0.0
    if row<=0:
        return (table[1,col]-table[0,col])/dcoord
    if row>=n-1:
        return (table[n-1,col]-table[n-2,col])/dcoord
    return 0.5*(table[row+1,col]-table[row-1,col])/dcoord

def antialias_weight(iw,nfft,dt,antialias,dtau):
    if antialias<=0.0 or iw<=0:
        return 1.0
    dtau=abs(dtau)*antialias
    if dtau<=1e-12 or not np.isfinite(dtau):
        return 1.0
    freq=iw/(nfft*dt)
    fnyq=0.5/dt
    fmax=0.5/dtau
    if fmax>=fnyq:
        return 1.0
    if fmax<=0.0:
        return 0.0
    fpass=0.8*fmax
    if freq<=fpass:
        return 1.0
    if freq>=fmax:
        return 0.0
    x=(freq-fpass)/(fmax-fpass)
    return 0.5*(1+math.cos(math.pi*x))

@dataclass
class Ctx:
    nt:int=512
    nh:int=256
    ns:int=256
    nsg:int=256
    nrg:int=256
    aper:int=64
    tap:int=10
    cmp:int=0
    nfft:int=2048
    nw:int=1025
    nsam:int=32
    dt:float=0.001
    h0:float=0.0
    dh:float=10.0
    s0:float=0.0
    ds:float=10.0
    sg0:float=0.0
    dsg:float=10.0
    rg0:float=0.0
    drg:float=10.0
    sdatum:float=120.0
    rdatum:float=120.0
    antialias:float=1.0
    stable:np.ndarray=None
    rtable:np.ndarray=None

def make_tables(ctx,vel=2200.0):
    xr=np.arange(ctx.nrg)*ctx.drg+ctx.rg0
    ctx.rtable=np.zeros((ctx.nrg,ctx.nrg),float)
    for cc in range(ctx.nrg):
        for c in range(ctx.nrg):
            dx=xr[cc]-xr[c]
            ctx.rtable[cc,c]=math.sqrt(ctx.rdatum*ctx.rdatum+dx*dx)/vel
    xs=np.arange(ctx.nsg)*ctx.dsg+ctx.sg0
    ctx.stable=np.zeros((ctx.nsg,ctx.nsg),float)
    for cc in range(ctx.nsg):
        for c in range(ctx.nsg):
            dx=xs[cc]-xs[c]
            ctx.stable[cc,c]=math.sqrt(ctx.sdatum*ctx.sdatum+dx*dx)/vel

class Kernel:
    def __init__(self,ctx):
        self.ctx=ctx
        self.cache={}
    def H(self,tau,iw):
        key=(iw,round(float(tau),12))
        if key not in self.cache:
            self.cache[key]=Hfast(float(tau),iw,self.ctx.nfft,self.ctx.dt,self.ctx.nsam)
        return self.cache[key]

    def receiver_tau(self,iw,is_,ih,ic):
        c=int((((self.ctx.s0+is_*self.ctx.ds) if self.ctx.cmp else 0.0)
               + self.ctx.h0 + ih*self.ctx.dh - self.ctx.rg0)/self.ctx.drg + 0.5)
        cc=int((((self.ctx.s0+is_*self.ctx.ds) if self.ctx.cmp else 0.0)
                + self.ctx.h0 + ic*self.ctx.dh - self.ctx.rg0)/self.ctx.drg + 0.5)
        if c<0 or c>=self.ctx.nrg or cc<0 or cc>=self.ctx.nrg:
            return None,None,None
        return self.ctx.rtable[cc,c],cc,c

    def receiver_full(self,iw,is_,ih,ic):
        ctx=self.ctx
        left=max(0,ih-ctx.aper); right=min(ctx.nh-1,ih+ctx.aper)
        if ic<left or ic>right:
            return 0j
        tau,cc,c=self.receiver_tau(iw,is_,ih,ic)
        if tau is None:
            return 0j
        coef=taper_weight(left,ic,right,ctx.tap,1)
        dist=ctx.rdatum*ctx.rdatum+(ic-ih)*(ic-ih)*ctx.dh*ctx.dh
        w=coef/math.pi*ctx.dh*ctx.rdatum*tau/dist
        slope=table_slope_first_index(ctx.rtable,cc,c,ctx.drg)
        aa=antialias_weight(iw,ctx.nfft,ctx.dt,ctx.antialias,slope*abs(ctx.dh))
        return w*aa*self.H(tau,iw)

    def receiver_H(self,iw,is_,ih,ic):
        tau,_,_=self.receiver_tau(iw,is_,ih,ic)
        return 0j if tau is None else self.H(tau,iw)

    def source0_tau(self,iw,ih,is_,ic):
        c=int((self.ctx.s0+is_*self.ctx.ds-self.ctx.sg0)/self.ctx.dsg + 0.5)
        cc=int((self.ctx.s0+ic*self.ctx.ds-self.ctx.sg0)/self.ctx.dsg + 0.5)
        if c<0 or c>=self.ctx.nsg or cc<0 or cc>=self.ctx.nsg:
            return None,None,None
        return self.ctx.stable[cc,c],cc,c

    def source0_full(self,iw,ih,is_,ic):
        ctx=self.ctx
        left=max(0,is_-ctx.aper); right=min(ctx.ns-1,is_+ctx.aper)
        if ic<left or ic>right:
            return 0j
        tau,cc,c=self.source0_tau(iw,ih,is_,ic)
        if tau is None:
            return 0j
        coef=taper_weight(left,ic,right,ctx.tap,1)
        dist=ctx.sdatum*ctx.sdatum+(ic-is_)*(ic-is_)*ctx.ds*ctx.ds
        w=coef/math.pi*ctx.ds*ctx.sdatum*tau/dist
        slope=table_slope_first_index(ctx.stable,cc,c,ctx.dsg)
        aa=antialias_weight(iw,ctx.nfft,ctx.dt,ctx.antialias,slope*abs(ctx.ds))
        return w*aa*self.H(tau,iw)

    def source0_H(self,iw,ih,is_,ic):
        tau,_,_=self.source0_tau(iw,ih,is_,ic)
        return 0j if tau is None else self.H(tau,iw)

    def source1_tau(self,iw,gm,q,p):
        ctx=self.ctx
        if q<0 or q>=gm["nloc"] or p<0 or p>=gm["nloc"]:
            return None,None,None,None,None
        is_=gm["isrc"][q]; ic=gm["isrc"][p]
        ih=gm["ih"][q]; hh=gm["ih"][p]
        if ih<0 or ih>=ctx.nh or hh<0 or hh>=ctx.nh:
            return None,None,None,None,None
        c=int((ctx.s0+is_*ctx.ds-ctx.sg0)/ctx.dsg + 0.5)
        cc=int((ctx.s0+ic*ctx.ds-ctx.sg0)/ctx.dsg + 0.5)
        if c<0 or c>=ctx.nsg or cc<0 or cc>=ctx.nsg:
            return None,None,None,None,None
        return ctx.stable[cc,c],cc,c,is_,ic

    def source1_full(self,iw,gm,q,p):
        ctx=self.ctx
        tau,cc,c,is_,ic=self.source1_tau(iw,gm,q,p)
        if tau is None:
            return 0j
        left=max(gm["sleft"],is_-gm["jump"]*ctx.aper)
        right=min(gm["sright"],is_+gm["jump"]*ctx.aper)
        if ic<left or ic>right:
            return 0j
        coef=taper_weight(left,ic,right,ctx.tap,gm["jump"])
        dist=ctx.sdatum*ctx.sdatum+(ic-is_)*(ic-is_)*ctx.ds*ctx.ds
        w=coef/math.pi*ctx.ds*ctx.sdatum*tau/dist
        slope=table_slope_first_index(ctx.stable,cc,c,ctx.dsg)
        aa=antialias_weight(iw,ctx.nfft,ctx.dt,ctx.antialias,slope*abs(gm["jump"]*ctx.ds))
        return w*aa*self.H(tau,iw)

    def source1_H(self,iw,gm,q,p):
        tau,_,_,_,_=self.source1_tau(iw,gm,q,p)
        return 0j if tau is None else self.H(tau,iw)

def build_cmp1_gather_map(ctx,ir):
    s=abs((ctx.ns-1)*ctx.ds)
    h=abs((ctx.nh-1)*ctx.dh)
    jump=1 if abs(ctx.ds)>=abs(ctx.dh) else int(ctx.dh/ctx.ds + 0.5)
    dr=abs(ctx.dh) if abs(ctx.ds)>=abs(ctx.dh) else abs(ctx.ds)
    r=ir*dr + (-1.0 if ctx.ds<=0 else 0.0)*s + (-1.0 if ctx.dh<=0 else 0.0)*h
    sleft=int((ir*dr + (-1.0 if ctx.ds<=0 else 0.0)*s + (0.0 if ctx.ds<=0 else -1.0)*h)/ctx.ds + 0.5)
    sright=int((ir*dr + (-1.0 if ctx.ds<=0 else 0.0)*s + (-1.0 if ctx.ds<=0 else 0.0)*h)/ctx.ds + 0.5)
    if sleft<0: sleft=0
    if sright>ctx.ns-1: sright=ctx.ns-1
    left0=int((r-sleft*ctx.ds)/ctx.dh + 0.5)
    if left0<0 or left0>ctx.nh-1: sleft+=1
    right0=int((r-sright*ctx.ds)/ctx.dh + 0.5)
    if right0<0 or right0>ctx.nh-1: sright-=1
    if sright<sleft:
        return dict(ir=ir,r=r,sleft=sleft,sright=sright,jump=jump,nloc=0,isrc=[],ih=[])
    isrc=list(range(sleft,sright+1,jump))
    ih=[int((r-is_*ctx.ds)/ctx.dh + 0.5) for is_ in isrc]
    return dict(ir=ir,r=r,sleft=sleft,sright=sright,jump=jump,nloc=len(isrc),isrc=isrc,ih=ih)

def direct_apply(nout,nin,full,x):
    y=np.zeros(nout,dtype=np.complex128)
    for i in range(nout):
        s=0j
        for j in range(nin):
            s += full(i,j)*x[j]
        y[i]=s
    return y

def bf_phase_residual(nout,nin,x,xin0,dxin,xout0,dxout,full,Hfunc,leaf_size=16,order=64):
    source=Tree1D(nin,xin0,dxin,leaf_size,order)
    target=Tree1D(nout,xout0,dxout,leaf_size,order)
    if source.depth!=target.depth:
        return direct_apply(nout,nin,full,x)
    L=source.depth
    mid=L//2
    sboxes=source.boxes
    tboxes=target.boxes

    sigma={}
    Aroot=target.root()
    A=tboxes[Aroot]
    out_center=A.center_index

    for bid in source.level_boxes(L):
        B=sboxes[bid]
        q=len(B.cheb_index)
        coeff=np.zeros(q,dtype=np.complex128)
        for j in range(B.begin,B.end):
            xcoord=xin0+j*dxin
            Lval=lagrange_all(B.cheb_coord,B.bary_w,xcoord)
            for t,in_node in enumerate(B.cheb_index):
                P=phase_ratio(Hfunc(out_center,j),Hfunc(out_center,in_node))
                coeff[t]+=Lval[t]*P*x[j]
        sigma[(Aroot,bid)]=coeff

    for level in range(0,mid):
        new={}
        for aid in target.level_boxes(level+1):
            A=tboxes[aid]
            apid=A.parent
            out_center=A.center_index
            for bpid in source.level_boxes(L-level-1):
                B=sboxes[bpid]
                q=len(B.cheb_index)
                out=np.zeros(q,dtype=np.complex128)
                for bcid in [B.left,B.right]:
                    if bcid<0: continue
                    old=sigma.get((apid,bcid))
                    if old is None: continue
                    Bc=sboxes[bcid]
                    for tp,in_child in enumerate(Bc.cheb_index):
                        xchild=Bc.cheb_coord[tp]
                        Lval=lagrange_all(B.cheb_coord,B.bary_w,xchild)
                        for t,in_parent in enumerate(B.cheb_index):
                            P=phase_ratio(Hfunc(out_center,in_child),Hfunc(out_center,in_parent))
                            out[t]+=Lval[t]*P*old[tp]
                new[(aid,bpid)]=out
        sigma=new

    gamma={}
    for (aid,bid),sig in sigma.items():
        A=tboxes[aid]; B=sboxes[bid]
        gam=np.zeros(len(A.cheb_index),dtype=np.complex128)
        for ta,iout in enumerate(A.cheb_index):
            s=0j
            for tb,iin in enumerate(B.cheb_index):
                s += full(iout,iin)*sig[tb]
            gam[ta]=s
        gamma[(aid,bid)]=gam

    for level in range(mid,L):
        new={}
        for aid in target.level_boxes(level+1):
            A=tboxes[aid]
            apid=A.parent
            Ap=tboxes[apid]
            for bpid in source.level_boxes(L-level-1):
                B=sboxes[bpid]
                out=np.zeros(len(A.cheb_index),dtype=np.complex128)
                for bcid in [B.left,B.right]:
                    if bcid<0: continue
                    old=gamma.get((apid,bcid))
                    if old is None: continue
                    Bc=sboxes[bcid]
                    in_center=Bc.center_index
                    for ta,out_child in enumerate(A.cheb_index):
                        rchild=A.cheb_coord[ta]
                        Lval=lagrange_all(Ap.cheb_coord,Ap.bary_w,rchild)
                        for tp,out_parent in enumerate(Ap.cheb_index):
                            P=phase_ratio(Hfunc(out_child,in_center),Hfunc(out_parent,in_center))
                            out[ta]+=Lval[tp]*P*old[tp]
                new[(aid,bpid)]=out
        gamma=new

    y=np.zeros(nout,dtype=np.complex128)
    source_root=source.root()
    Broot=sboxes[source_root]
    in_center=Broot.center_index

    for aid in target.level_boxes(L):
        A=tboxes[aid]
        gam=gamma.get((aid,source_root))
        if gam is None:
            continue
        for i in range(A.begin,A.end):
            rcoord=xout0+i*dxout
            Lval=lagrange_all(A.cheb_coord,A.bary_w,rcoord)
            s=0j
            for t,out_node in enumerate(A.cheb_index):
                P=phase_ratio(Hfunc(i,in_center),Hfunc(out_node,in_center))
                s += Lval[t]*P*gam[t]
            y[i]=s
    return y


def run_one_case(iw, order, leaf_size, n=256, aper=64, tap=10, antialias=1.0):
    """
    Validate three complete current extrapolation kernels with order < leaf_size.
    Returns relative L2 errors for receiver, source cmp=0, and source cmp=1 local gather.
    """
    if not (order < leaf_size):
        raise ValueError(f"Expected order < leaf_size, got order={order}, leaf_size={leaf_size}")

    ctx = Ctx(nh=n, ns=n, nrg=n, nsg=n, nfft=2048, nw=1025,
              aper=aper, tap=tap, antialias=antialias)
    make_tables(ctx, vel=2200.0)
    ker = Kernel(ctx)

    np.random.seed(2026)

    # 1. receiver-side
    is_fixed = min(10, ctx.ns - 1)
    xrec = np.exp(1j * 0.05 * np.arange(ctx.nh)) * (1.0 + 0.1 * np.sin(0.07 * np.arange(ctx.nh)))
    full_rec = lambda ih, ic: ker.receiver_full(iw, is_fixed, ih, ic)
    H_rec = lambda ih, ic: ker.receiver_H(iw, is_fixed, ih, ic)
    yd = direct_apply(ctx.nh, ctx.nh, full_rec, xrec)
    yb = bf_phase_residual(ctx.nh, ctx.nh, xrec,
                           ctx.h0, ctx.dh, ctx.h0, ctx.dh,
                           full_rec, H_rec, leaf_size, order)
    e_rec = relerr(yb, yd)

    # 2. source-side cmp=0
    ih_fixed = min(20, ctx.nh - 1)
    xsrc = np.exp(1j * 0.04 * np.arange(ctx.ns)) * (1.0 + 0.2 * np.cos(0.05 * np.arange(ctx.ns)))
    full_s0 = lambda is_, ic: ker.source0_full(iw, ih_fixed, is_, ic)
    H_s0 = lambda is_, ic: ker.source0_H(iw, ih_fixed, is_, ic)
    yd0 = direct_apply(ctx.ns, ctx.ns, full_s0, xsrc)
    yb0 = bf_phase_residual(ctx.ns, ctx.ns, xsrc,
                            ctx.s0, ctx.ds, ctx.s0, ctx.ds,
                            full_s0, H_s0, leaf_size, order)
    e_s0 = relerr(yb0, yd0)

    # 3. source-side cmp=1 local common-receiver gather
    ctx1 = Ctx(nh=n, ns=n, nrg=n, nsg=n, nfft=2048, nw=1025,
               aper=aper, tap=tap, antialias=antialias, cmp=1)
    make_tables(ctx1, vel=2200.0)
    ker1 = Kernel(ctx1)

    s = abs((ctx1.ns - 1) * ctx1.ds)
    h = abs((ctx1.nh - 1) * ctx1.dh)
    dr = abs(ctx1.dh) if abs(ctx1.ds) >= abs(ctx1.dh) else abs(ctx1.ds)
    nr = int((s + h) / dr + 1.5)
    ir = nr // 2
    gm = build_cmp1_gather_map(ctx1, ir)

    x1 = np.exp(1j * 0.04 * np.arange(gm["nloc"])) * (1.0 + 0.2 * np.cos(0.05 * np.arange(gm["nloc"])))
    full_1 = lambda q, p: ker1.source1_full(iw, gm, q, p)
    H_1 = lambda q, p: ker1.source1_H(iw, gm, q, p)

    yd1 = direct_apply(gm["nloc"], gm["nloc"], full_1, x1)
    yb1 = bf_phase_residual(gm["nloc"], gm["nloc"], x1,
                            0.0, 1.0, 0.0, 1.0,
                            full_1, H_1, leaf_size, order)
    e_s1 = relerr(yb1, yd1)

    return {
        "iw": iw,
        "order": order,
        "leaf_size": leaf_size,
        "receiver": e_rec,
        "source_cmp0": e_s0,
        "source_cmp1": e_s1,
        "cmp1_ir": ir,
        "cmp1_nloc": gm["nloc"],
    }


def run_sweep():
    """
    Sweep only valid configurations with order < leaf_size.
    This intentionally avoids invalid cases such as order=64, leaf_size=16.
    """
    configs = [
        (8, 16),
        (16, 32),
        (32, 64),
        (60, 64),
    ]
    frequencies = [1, 5, 10, 20, 40, 80]

    results = []
    for order, leaf in configs:
        for iw in frequencies:
            res = run_one_case(iw=iw, order=order, leaf_size=leaf, n=128)
            results.append(res)
            print(
                f"order={order:2d} leaf={leaf:2d} iw={iw:3d} | "
                f"rec={res['receiver']:.6e} "
                f"src0={res['source_cmp0']:.6e} "
                f"src1={res['source_cmp1']:.6e}",
                flush=True
            )

    return results


def write_csv(results, path):
    with open(path, "w", encoding="utf-8") as f:
        f.write("order,leaf_size,iw,receiver,source_cmp0,source_cmp1,cmp1_ir,cmp1_nloc\n")
        for r in results:
            f.write(
                f"{r['order']},{r['leaf_size']},{r['iw']},"
                f"{r['receiver']:.12e},{r['source_cmp0']:.12e},{r['source_cmp1']:.12e},"
                f"{r['cmp1_ir']},{r['cmp1_nloc']}\n"
            )


if __name__ == "__main__":
    results = run_sweep()
    csv_path = "/mnt/data/order_lt_leaf_current_kernel_validation_results.csv"
    write_csv(results, csv_path)

    print("")
    print(f"Saved CSV results to: {csv_path}")
    print("")
    print("Interpretation:")
    print("- These tests all satisfy order < leaf_size.")
    print("- Low frequencies can reach small error.")
    print("- High frequencies show large errors for the current complete K = W * AA * H_fast kernel.")
    print("- Therefore the current phase-residual butterfly guide must not claim full-kernel validation at all frequencies.")

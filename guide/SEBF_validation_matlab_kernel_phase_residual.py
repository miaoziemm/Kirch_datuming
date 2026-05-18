import numpy as np

def bary_weights(nodes):
    nodes=np.asarray(nodes,float)
    q=len(nodes)
    w=np.ones(q)
    for j in range(q):
        for k in range(q):
            if k!=j:
                w[j]/=(nodes[j]-nodes[k])
    return w

def lagrange_all(x,nodes,w):
    nodes=np.asarray(nodes,float); w=np.asarray(w,float)
    d=x-nodes
    k=np.where(np.abs(d)<1e-13)[0]
    if k.size:
        L=np.zeros_like(nodes); L[k[0]]=1.0; return L
    tmp=w/d
    return tmp/tmp.sum()

def build_Lmat(points,nodes,w):
    return np.vstack([lagrange_all(p,nodes,w) for p in points])

class Box:
    pass

def build_tree(grid, leaf_n, cheb_std):
    n=len(grid)
    L=int(round(np.log2(n/leaf_n)))
    assert abs(2**L*leaf_n-n)<1e-12
    T=[]
    for lev in range(L+1):
        nbox=2**lev; box_n=n//nbox
        boxes=[]
        for b in range(nbox):
            i1=b*box_n; i2=(b+1)*box_n-1
            xL=grid[i1]; xR=grid[i2]
            box=Box(); box.i1=i1; box.i2=i2
            box.x0=0.5*(xL+xR); box.w=xR-xL
            box.cheb=box.x0+0.5*box.w*cheb_std
            boxes.append(box)
        T.append(boxes)
    return T

# continuous phase/amplitude equivalent to MATLAB example
def tau_func(s,r,v=3000.0,dz=100.0):
    return np.sqrt((s-r)**2+dz**2)/v

def amp_func(s,r,dz=100.0):
    rr=np.sqrt((s-r)**2+dz**2)
    return 1/np.sqrt(rr)*dz/rr

def validate(n=256, leaf_n=16, p=20, dx=5.0, omega=15*2*np.pi):
    v=3000.0; dz=100.0
    grid=np.arange(1,n+1)*dx
    f=np.zeros(n,dtype=np.complex128)
    f[n//2-1]=1; f[n//2]=1
    f=f*1j*omega/v
    # direct
    S,R=np.meshgrid(grid,grid,indexing='ij')
    tt=tau_func(S,R,v,dz)
    amp=amp_func(S,R,dz)
    K=amp*np.exp(1j*omega*tt)
    direct=(K*f[:,None]).sum(axis=0)
    k=np.arange(1,p+1)
    cheb_std=np.cos((2*k-1)*np.pi/(2*p))
    TR=build_tree(grid,leaf_n,cheb_std); TS=build_tree(grid,leaf_n,cheb_std)
    L=int(round(np.log2(n/leaf_n)))
    # init sigma shape (nbox_s at leaf, nbox_r root=1, p)
    nbox_s=len(TS[L]); nbox_r=len(TR[0])
    sigma=np.zeros((nbox_s,nbox_r,p),dtype=np.complex128)
    for ibs,B in enumerate(TS[L]):
        nodes=B.cheb; w=bary_weights(nodes)
        Brange=grid[B.i1:B.i2+1]
        Lg=build_Lmat(Brange,nodes,w)
        for ibr,A in enumerate(TR[0]):
            r0=A.x0
            phase1=tau_func(Brange,r0,v,dz)
            temp=Lg.T@(np.exp(1j*omega*phase1)*f[B.i1:B.i2+1])
            phase2=tau_func(nodes,r0,v,dz)
            sigma[ibs,ibr,:]=temp*np.exp(-1j*omega*phase2)
    lv_max=L//2
    sigma_prev=sigma
    for lv in range(1,lv_max+1):
        boxes_r_now=TR[lv]; boxes_s_now=TS[L-lv]
        boxes_r_prev=TR[lv-1]; boxes_s_prev=TS[L-lv+1]
        sigma_curr=np.zeros((len(boxes_s_now),len(boxes_r_now),p),dtype=np.complex128)
        for ibr,A in enumerate(boxes_r_now):
            r0=A.x0; ibr_parent=(ibr)//2
            for ibs,B in enumerate(boxes_s_now):
                nodesB=B.cheb; wB=bary_weights(nodesB)
                acc=np.zeros(p,dtype=np.complex128)
                for child in [2*ibs,2*ibs+1]:
                    Bc=boxes_s_prev[child]
                    child_nodes=Bc.cheb
                    Lmat=build_Lmat(child_nodes,nodesB,wB)
                    sigma_child=sigma_prev[child,ibr_parent,:]
                    phase_child=tau_func(child_nodes,r0,v,dz)
                    tmp=np.exp(1j*omega*phase_child)*sigma_child
                    acc += Lmat.T @ tmp
                phase_B=tau_func(nodesB,r0,v,dz)
                sigma_curr[ibs,ibr,:]=acc*np.exp(-1j*omega*phase_B)
        sigma_prev=sigma_curr
    sigma=sigma_prev
    # switch
    mid_lv=lv_max
    boxes_r_mid=TR[mid_lv]; boxes_s_mid=TS[L-mid_lv]
    sigma_sw=np.zeros((len(boxes_s_mid),len(boxes_r_mid),p),dtype=np.complex128)
    for ibr,A in enumerate(boxes_r_mid):
        r_nodes=A.cheb
        for ibs,B in enumerate(boxes_s_mid):
            s_nodes=B.cheb
            sigma_old=sigma[ibs,ibr,:]
            S2,R2=np.meshgrid(s_nodes,r_nodes,indexing='ij') # shape p,p? rows source, cols target
            # want K_switch(target it, source jt)
            Ksw=(amp_func(S2.T,R2.T,dz)*np.exp(1j*omega*tau_func(S2.T,R2.T,v,dz)))
            sigma_sw[ibs,ibr,:]=Ksw @ sigma_old
    sigma_prev=sigma_sw
    # second half
    for lv in range(mid_lv+1,L+1):
        boxes_r_now=TR[lv]; boxes_s_now=TS[L-lv]
        boxes_r_prev=TR[lv-1]; boxes_s_prev=TS[L-lv+1]
        sigma_curr=np.zeros((len(boxes_s_now),len(boxes_r_now),p),dtype=np.complex128)
        for ibr,A in enumerate(boxes_r_now):
            r_nodes=A.cheb; ibr_parent=ibr//2; Ap=boxes_r_prev[ibr_parent]; rp_nodes=Ap.cheb
            wAp=bary_weights(rp_nodes); Lmat=build_Lmat(r_nodes,rp_nodes,wAp)
            for ibs,B in enumerate(boxes_s_now):
                acc=np.zeros(p,dtype=np.complex128)
                for child in [2*ibs,2*ibs+1]:
                    Bc=boxes_s_prev[child]; s0Bc=Bc.x0
                    sigma_child=sigma_prev[child,ibr_parent,:]
                    fac_r=np.exp(1j*omega*tau_func(s0Bc,r_nodes,v,dz))
                    fac_rp=np.exp(-1j*omega*tau_func(s0Bc,rp_nodes,v,dz))
                    tmp=fac_rp*sigma_child
                    acc += fac_r*(Lmat @ tmp)
                sigma_curr[ibs,ibr,:]=acc
        sigma_prev=sigma_curr
    sigma=sigma_prev
    # termination
    boxes_r_leaf=TR[L]; Broot=TS[0][0]; s0B=Broot.x0
    u=np.zeros(n,dtype=np.complex128)
    for ibr,A in enumerate(boxes_r_leaf):
        Arange=grid[A.i1:A.i2+1]
        r_nodes=A.cheb; wA=bary_weights(r_nodes)
        LgA=build_Lmat(Arange,r_nodes,wA)
        sigma_leaf=sigma[0,ibr,:]
        fac_rt=np.exp(-1j*omega*tau_func(s0B,r_nodes,v,dz))
        fac_r=np.exp(1j*omega*tau_func(s0B,Arange,v,dz))
        uA=fac_r*(LgA@(fac_rt*sigma_leaf))
        u[A.i1:A.i2+1]=uA
    rel=np.linalg.norm(u-direct)/np.linalg.norm(direct)
    return rel,u,direct

if __name__=='__main__':
    for p in [8,12,16,20,24,32]:
        rel,_,_=validate(n=256,leaf_n=16,p=p)
        print('p',p,'rel',rel)

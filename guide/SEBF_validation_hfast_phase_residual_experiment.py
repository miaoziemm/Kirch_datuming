import numpy as np
from validate_phase_residual_bf import bary_weights, lagrange_all, build_Lmat, build_tree, amp_func, tau_func

def Hfast(tau, iw=60, nfft=1024, dt=0.001, nsam=32):
    tau=np.asarray(tau,dtype=float)
    # vectorized output complex
    G=np.zeros(tau.shape+(nsam,),dtype=float)
    for k in range(1,nsam):
        x=(tau+k*dt)/tau
        G[...,k]=np.sqrt(x*x-1.0)
    D=G[...,1:]-G[...,:-1] # 0..nsam-2 diff first
    # second derivative from nsam-2 down to 1: F[k]=D[k]-D[k-1], F[0]=D[0]
    F=D.copy()
    F[...,1:]=D[...,1:]-D[...,:-1]
    q=tau/dt; m=np.ceil(q).astype(int); delta=m-q
    omega=2*np.pi*iw/nfft
    # sumF
    karr=np.arange(nsam-1)
    expk=np.exp(-1j*omega*karr)
    sumF=(F/dt*expk).sum(axis=-1)
    interp=(1-delta)*np.exp(-1j*omega*m)+delta*np.exp(-1j*omega*(m-1))
    return interp*sumF

def phase_factor(Ha,Hb):
    z=Ha*np.conj(Hb); mag=np.abs(z)
    return np.where(mag<1e-20,1+0j,z/mag)

def validate(n=256, leaf_n=16, p=20, dx=5.0, iw=60):
    v=3000.0; dz=100.0
    grid=np.arange(1,n+1)*dx
    f=np.zeros(n,dtype=np.complex128); f[n//2-1]=1; f[n//2]=1
    # direct current-like kernel K = amp * Hfast(tau)
    S,R=np.meshgrid(grid,grid,indexing='ij')
    tau=tau_func(S,R,v,dz); amp=amp_func(S,R,dz)
    K=amp*Hfast(tau,iw=iw)
    direct=(K*f[:,None]).sum(axis=0)
    k=np.arange(1,p+1); cheb_std=np.cos((2*k-1)*np.pi/(2*p))
    TR=build_tree(grid,leaf_n,cheb_std); TS=build_tree(grid,leaf_n,cheb_std)
    L=int(round(np.log2(n/leaf_n)))
    # init with phase H ratios, not full amp
    sigma=np.zeros((len(TS[L]),len(TR[0]),p),dtype=np.complex128)
    for ibs,B in enumerate(TS[L]):
        nodes=B.cheb; w=bary_weights(nodes); Brange=grid[B.i1:B.i2+1]
        Lg=build_Lmat(Brange,nodes,w)
        for ibr,A in enumerate(TR[0]):
            r0=A.x0
            H1=Hfast(tau_func(Brange,r0,v,dz),iw=iw)
            H2=Hfast(tau_func(nodes,r0,v,dz),iw=iw)
            # coeff[t] = sum_j L[j,t] P((r0,sj),(r0,st)) xj
            for t in range(p):
                P=phase_factor(H1,H2[t])
                sigma[ibs,ibr,t]=np.sum(Lg[:,t]*P*f[B.i1:B.i2+1])
    lv_max=L//2; sigma_prev=sigma
    for lv in range(1,lv_max+1):
        boxes_r_now=TR[lv]; boxes_s_now=TS[L-lv]
        boxes_r_prev=TR[lv-1]; boxes_s_prev=TS[L-lv+1]
        sigma_curr=np.zeros((len(boxes_s_now),len(boxes_r_now),p),dtype=np.complex128)
        for ibr,A in enumerate(boxes_r_now):
            r0=A.x0; ibr_parent=ibr//2
            for ibs,B in enumerate(boxes_s_now):
                nodesB=B.cheb; wB=bary_weights(nodesB); acc=np.zeros(p,dtype=np.complex128)
                H_parent=Hfast(tau_func(nodesB,r0,v,dz),iw=iw)
                L_targets_cache=[]
                for child in [2*ibs,2*ibs+1]:
                    Bc=boxes_s_prev[child]; child_nodes=Bc.cheb
                    Lmat=build_Lmat(child_nodes,nodesB,wB) # child p x parent p
                    sigma_child=sigma_prev[child,ibr_parent,:]
                    H_child=Hfast(tau_func(child_nodes,r0,v,dz),iw=iw)
                    for t in range(p):
                        P=phase_factor(H_child,H_parent[t]) # vector over t'
                        acc[t]+=np.sum(Lmat[:,t]*P*sigma_child)
                sigma_curr[ibs,ibr,:]=acc
        sigma_prev=sigma_curr
    sigma=sigma_prev
    # switch full K
    mid=lv_max; sigma_sw=np.zeros_like(sigma)
    for ibr,A in enumerate(TR[mid]):
        r_nodes=A.cheb
        for ibs,B in enumerate(TS[L-mid]):
            s_nodes=B.cheb; sigma_old=sigma[ibs,ibr,:]
            S2,R2=np.meshgrid(s_nodes,r_nodes,indexing='xy') # R rows?, easier
            # R2 shape p,p with rows r? Actually meshgrid x=s, y=r gives shape (p,p), rows r, cols s
            tau2=tau_func(S2,R2,v,dz); amp2=amp_func(S2,R2,dz)
            Ksw=amp2*Hfast(tau2,iw=iw)
            sigma_sw[ibs,ibr,:]=Ksw@sigma_old
    sigma_prev=sigma_sw
    # second half
    for lv in range(mid+1,L+1):
        boxes_r_now=TR[lv]; boxes_s_now=TS[L-lv]
        boxes_r_prev=TR[lv-1]; boxes_s_prev=TS[L-lv+1]
        sigma_curr=np.zeros((len(boxes_s_now),len(boxes_r_now),p),dtype=np.complex128)
        for ibr,A in enumerate(boxes_r_now):
            r_nodes=A.cheb; ibr_parent=ibr//2; Ap=boxes_r_prev[ibr_parent]; rp_nodes=Ap.cheb
            wAp=bary_weights(rp_nodes); Lmat=build_Lmat(r_nodes,rp_nodes,wAp)
            for ibs,B in enumerate(boxes_s_now):
                acc=np.zeros(p,dtype=np.complex128)
                for child in [2*ibs,2*ibs+1]:
                    Bc=boxes_s_prev[child]; s0=Bc.x0; sigma_child=sigma_prev[child,ibr_parent,:]
                    H_r=Hfast(tau_func(s0,r_nodes,v,dz),iw=iw) # len p
                    H_rp=Hfast(tau_func(s0,rp_nodes,v,dz),iw=iw)
                    for t in range(p):
                        P=phase_factor(H_r[t],H_rp) # vector over t'
                        acc[t]+=np.sum(Lmat[t,:]*P*sigma_child)
                sigma_curr[ibs,ibr,:]=acc
        sigma_prev=sigma_curr
    sigma=sigma_prev
    # termination
    u=np.zeros(n,dtype=np.complex128); Broot=TS[0][0]; s0=Broot.x0
    for ibr,A in enumerate(TR[L]):
        Arange=grid[A.i1:A.i2+1]; r_nodes=A.cheb; wA=bary_weights(r_nodes)
        Lg=build_Lmat(Arange,r_nodes,wA); sig=sigma[0,ibr,:]
        H_real=Hfast(tau_func(s0,Arange,v,dz),iw=iw)
        H_node=Hfast(tau_func(s0,r_nodes,v,dz),iw=iw)
        uA=np.zeros(len(Arange),dtype=np.complex128)
        for i in range(len(Arange)):
            P=phase_factor(H_real[i],H_node)
            uA[i]=np.sum(Lg[i,:]*P*sig)
        u[A.i1:A.i2+1]=uA
    return np.linalg.norm(u-direct)/np.linalg.norm(direct),u,direct

if __name__=='__main__':
    for p in [8,12,16,20,24,32,40]:
        rel,_,_=validate(n=256,leaf_n=16,p=p,iw=60)
        print('p',p,'rel',rel)

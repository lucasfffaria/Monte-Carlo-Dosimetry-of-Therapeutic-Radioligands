"""Phantom data (ICRP 110 AM/AF, standing and seated) for MCNP6 lattice models."""
import numpy as np, json, os
import paths, phantom_cache
BASE=paths.ICRP110
VOX={'AM':(0.2137,0.2137,0.8),'AF':(0.1775,0.1775,0.484)}
def organ_table(sex):
    d={}
    for line in open(BASE/sex/f"{sex}_organs.dat",encoding='latin-1'):
        p=line.split()
        if len(p)>3 and p[0].isdigit():
            d[int(p[0])]=dict(name=' '.join(p[1:-2]),medium=int(p[-2]),rho=float(p[-1]))
    d[141]=dict(name='Skin (top/bottom slice)',medium=27,rho=1.09)
    return d
def media_table(sex):
    lines=open(BASE/sex/f"{sex}_media.dat",encoding='latin-1').read().splitlines()
    Z=[1,6,7,8,11,12,15,16,17,19,20,26,53]
    med={}
    for line in lines:
        p=line.split()
        if len(p)>13 and p[0].isdigit():
            fr=[float(x) for x in p[-13:]]
            med[int(p[0])]=dict(name=' '.join(p[1:-13]),frac=dict(zip(Z,fr)))
    return med
class Phantom:
    """kind: 'std' or 'seat'; arr indexed [k,j,i] (slice,row,col); lattice index of arr[0,0,0] = (i0,j0,k0)"""
    def __init__(self,sex,kind):
        self.sex=sex; self.kind=kind
        self.dx,self.dy,self.dz=VOX[sex]
        if kind=='std':
            self.arr=phantom_cache.std_array(sex)
            self.i0=self.j0=self.k0=0
            self.meta={}
        else:
            self.arr=phantom_cache.seated_array(sex).astype(np.int32)
            self.meta=phantom_cache.seated_info(sex)
            self.i0,self.j0,self.k0=self.meta['origin_idx']
        self.nk,self.nj,self.ni=self.arr.shape
        self.organs=organ_table(sex); self.media=media_table(sex)
        ids,cnt=np.unique(self.arr,return_counts=True)
        self.count={int(i):int(c) for i,c in zip(ids,cnt) if i>0}
        self.vvox=self.dx*self.dy*self.dz
    # local (lattice-frame) coordinates of the lattice element (i,j,k): [i*dx,(i+1)*dx]...
    def extent(self):
        return (self.i0*self.dx,(self.i0+self.ni)*self.dx,self.j0*self.dy,(self.j0+self.nj)*self.dy,self.k0*self.dz,(self.k0+self.nk)*self.dz)
    def body_box(self,mask=None):
        m=self.arr>0 if mask is None else mask
        k,j,i=np.nonzero(m)
        return ((i.min()+self.i0)*self.dx,(i.max()+1+self.i0)*self.dx,(j.min()+self.j0)*self.dy,(j.max()+1+self.j0)*self.dy,(k.min()+self.k0)*self.dz,(k.max()+1+self.k0)*self.dz)
    def mass(self,oid):
        return self.count.get(oid,0)*self.vvox*self.organs[oid]['rho']
    def volume(self,oid):
        return self.count.get(oid,0)*self.vvox
    def write_fill(self,fn,cellno,mat,rho,u,uoff,surf6,header_comment=''):
        """write MCNP lattice cell card (LAT=1) with fully specified FILL; air voxels -> own universe u"""
        a=self.arr
        tok=[]
        flat=a.reshape(-1)          # i fastest, then j, then k  == MCNP FILL order
        univ=np.where(flat>0,flat+uoff,u).astype(np.int64)
        # run-length encode
        change=np.flatnonzero(np.diff(univ))+1
        starts=np.concatenate([[0],change]); ends=np.concatenate([change,[univ.size]])
        out=[]
        out.append(f"c {header_comment}")
        s=' '.join(str(x) for x in surf6)
        out.append(f"{cellno} {mat} -{rho} {s} lat=1 u={u} imp:p=1")
        out.append(f"     fill={self.i0}:{self.i0+self.ni-1} {self.j0}:{self.j0+self.nj-1} {self.k0}:{self.k0+self.nk-1}")
        line='     '
        for st,en in zip(starts,ends):
            v=univ[st]; n=en-st
            t=f"{v}" if n==1 else (f"{v} {v}" if n==2 else f"{v} {n-1}r")
            if len(line)+len(t)+1>78:
                out.append(line.rstrip()); line='     '
            line+=t+' '
        if line.strip(): out.append(line.rstrip())
        with open(fn,'w',newline='\r\n') as f:
            f.write('\n'.join(out)+'\n')
        return len(starts)

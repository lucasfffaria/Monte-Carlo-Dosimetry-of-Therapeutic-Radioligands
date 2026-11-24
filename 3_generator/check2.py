"""Extra checks for v2 decks (slab lattices):
 1. every filled container lies strictly inside the FILL index range of its lattice
 2. source emulation: for every SDEF entry, map its sampling box into lattice index space and check
    that it covers all voxels of the organ in that slab; MC estimate of acceptance vs expected;
    per-organ uniformity of probability per voxel across slabs
 3. tally chains: organ cells exist in the lattices filling the listed containers"""
import sys, os, re, numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import validate
def tr_of(c):
    t=c['tr']; o=np.array(t[:3])
    B=np.eye(3) if len(t)==3 else np.array(t[3:12]).reshape(3,3)
    return o,B
def lattice_info(cells,surfs,lat):
    c=cells[lat]; g=[int(x) for x in c['geo']]
    # geo: -s1 s2 -s3 s4 -s5 s6 with s1: px dx, s2: px 0 ...
    d=[surfs[abs(g[0])][1][0],surfs[abs(g[2])][1][0],surfs[abs(g[4])][1][0]]
    lo=[surfs[abs(g[1])][1][0],surfs[abs(g[3])][1][0],surfs[abs(g[5])][1][0]]
    assert lo==[0,0,0]
    i1,i2,j1,j2,k1,k2,vals=c['fillarr']
    arr=vals.reshape(k2-k1+1,j2-j1+1,i2-i1+1)
    return np.array(d),(i1,j1,k1),arr
def to_univ(p,o,B):
    # MCNP: x_cell = B x_univ + o  (B orthogonal) -> x_univ = B^T (x_cell - o)
    return (np.asarray(p)-o)@B.T
def check(fn,nmc0=4000,verbose=False):
    title,cells,surfs,dc=validate.parse(fn)
    msgs=[]
    # ---- 1. containers vs lattice ranges
    conts=[n for n,c in cells.items() if 'fill' in c and 'u' not in c]
    latinfo={}
    for n in conts:
        c=cells[n]; lat=[m for m,x in cells.items() if x.get('u')==c['fill'] and 'fillarr' in x]
        if not lat: continue
        lat=lat[0]
        if lat not in latinfo: latinfo[lat]=lattice_info(cells,surfs,lat)
        d,(i1,j1,k1),arr=latinfo[lat]
        assert len(c['geo'])==1 and c['geo'][0].startswith('-'),('container geo',n)
        b=surfs[int(c['geo'][0][1:])][1]
        o,B=tr_of(c)
        corners=np.array([[b[i],b[2+j],b[4+k]] for i in (0,1) for j in (0,1) for k in (0,1)])
        u=to_univ(corners,o,B); lo=u.min(0); hi=u.max(0)
        idx_lo=lo/d; idx_hi=hi/d
        rng_lo=np.array([i1,j1,k1]); rng_hi=rng_lo+np.array(arr.shape[::-1])
        margin=min((idx_lo-rng_lo).min(),(rng_hi-idx_hi).min())
        if margin<=0: msgs.append(f'container {n} leaves FILL range of lattice {lat} (margin {margin:.2e} elements)')
        c['_lat']=lat
    # ---- 2. source
    dists={}
    for cd in dc:
        t=cd.split(); m=re.match(r'^(SI|SP|DS|SB)(\d+)$',t[0].upper())
        if m: dists.setdefault(int(m.group(2)),{})[m.group(1)]=t[1:]
    sd=[x for x in dc if x.upper().startswith('SDEF')]
    stat=None
    if sd and 'CEL=D1' in sd[0].upper():
        paths=re.findall(r'\(([^)]*)\)',' '.join(dists[1]['SI'][1:]))
        p=np.array([float(x) for x in dists[1]['SP']])
        ds={k:[int(x) for x in dists[k]['DS'][1:]] for k in (2,3,4,5)}
        def dval(n):
            si=dists[n]['SI']
            if si[0].upper()=='H': return ('H',float(si[1]),float(si[2]))
            return ('L',float(si[1]))
        rng=np.random.default_rng(7)
        per_organ={}; effs=[]; bad=0
        for e,path in enumerate(paths):
            parts=[x.strip() for x in path.split('<')]
            cont=int(parts[-1].lstrip('-')); lat=int(parts[-2])
            c=cells[cont]; assert c['_lat']==lat,('path lattice/container mismatch',path)
            d,(i1,j1,k1),arr=latinfo[lat]; o,B=tr_of(c)
            cb=surfs[int(c['geo'][0][1:])][1]
            X=[dval(ds[k][e]) for k in (3,4,5)]
            if not parts[-1].startswith('-'):
                # v4 element sampling: PDS = organ level, X/Y/Z must span exactly one element
                assert len(parts)==3,('element-sampled path must start at the organ',path)
                okx=all(x[0]=='H' and abs(x[1])<1e-9 and abs(x[2]-dd)<1e-4 for x,dd in zip(X,d))
                if not okx: msgs.append(f'source entry {e} {path}: X/Y/Z not the element box {X}'); bad+=1
                oc=int(parts[0]); uo=cells[oc]['u']; ntot=int((arr==uo).sum())
                if ntot==0: msgs.append(f'source entry {e} {path}: organ not in lattice'); bad+=1
                per_organ.setdefault(oc,[]).append((p[e],ntot,ds[2][e]))
                effs.append((ntot/arr.size,None,path,0))
                continue
            assert parts[-1].startswith('-')
            # MC acceptance (sample size adapted to the expected efficiency)
            nmc=nmc0
            if len(parts)==3 and all(x[0]=='H' for x in X):
                vol=np.prod([x[2]-x[1] for x in X]); ntot0=int((arr==cells[int(parts[0])]['u']).sum())
                nmc=int(min(400000,max(nmc0,60/(ntot0*np.prod(d)/vol))))
            pts=np.empty((nmc,3))
            for a,x in enumerate(X):
                pts[:,a]=x[1] if x[0]=='L' else rng.uniform(x[1],x[2],nmc)
            inc=(pts[:,0]>cb[0])&(pts[:,0]<cb[1])&(pts[:,1]>cb[2])&(pts[:,1]<cb[3])&(pts[:,2]>cb[4])&(pts[:,2]<cb[5])
            u=to_univ(pts,o,B); ix=np.floor(u/d).astype(int)-np.array([i1,j1,k1])
            okr=(ix>=0).all(1)&(ix[:,0]<arr.shape[2])&(ix[:,1]<arr.shape[1])&(ix[:,2]<arr.shape[0])
            assert (okr|~inc).all(),('point in container outside lattice range',path)
            if len(parts)==3:
                oc=int(parts[0]); uo=cells[oc]['u']
                acc=np.zeros(nmc,bool); sel=inc&okr
                acc[sel]=arr[ix[sel,2],ix[sel,1],ix[sel,0]]==uo
                # exact coverage: organ voxels of this lattice inside the box (index space)
                if all(x[0]=='H' for x in X):
                    corners=np.array([[X[0][1+i],X[1][1+j],X[2][1+k]] for i in (0,1) for j in (0,1) for k in (0,1)])
                    uu=to_univ(corners,o,B); lo=np.round(uu.min(0)/d,2)-[i1,j1,k1]; hi=np.round(uu.max(0)/d,2)-[i1,j1,k1]
                    lo=np.maximum(np.floor(lo).astype(int),0); hi=np.ceil(hi).astype(int)
                    sub=arr[lo[2]:hi[2],lo[1]:hi[1],lo[0]:hi[0]]
                    ninbox=int((sub==uo).sum()); ntot=int((arr==uo).sum())
                    if ninbox!=ntot: msgs.append(f'source entry {e} {path}: box covers {ninbox} of {ntot} organ voxels'); bad+=1
                    per_organ.setdefault(oc,[]).append((p[e],ntot,ds[2][e]))
                    eff_exp=ntot/np.prod(hi-lo)
                else: eff_exp=None
                eff=acc.mean(); effs.append((eff,eff_exp,path,nmc))
                if eff==0: msgs.append(f'source entry {e} {path}: zero acceptance in {nmc} trials'); bad+=1
            else:
                acc=inc&okr
                if acc.mean()<0.99: msgs.append(f'source entry {e} {path}: point/line acceptance {acc.mean():.3f}'); bad+=1
        # uniformity per organ: p/ntot constant across slabs
        worst=0
        for oc,l in per_organ.items():
            r=np.array([a/b for a,b,_ in l]); worst=max(worst,r.max()/r.min()-1)
        # region probabilities
        reg={}
        for e in range(len(paths)): reg[ds[2][e]]=reg.get(ds[2][e],0)+p[e]
        regdev=max(abs(v*len(reg)-1) for v in reg.values())
        if worst>1e-4: msgs.append(f'non-uniform probability per voxel across slabs: {worst:.2e}')
        if regdev>1e-4: msgs.append(f'region probabilities not equal: {regdev:.2e}')
        e_arr=np.array([x[0] for x in effs])
        stat=dict(entries=len(paths),min_mc_eff=float(e_arr.min()) if len(e_arr) else None,bad=bad,
                  uniform_dev=worst,region_dev=regdev)
        # MC vs exact efficiency for a few entries
        dev=[(abs(a-b)/np.sqrt(max(b*(1-b)/n,1.0/n**2)),pth) for a,b,pth,n in effs if b]
        stat['worst_sigma']=max(dev)[0] if dev else None; stat['worst_sigma_entry']=max(dev)[1] if dev else None
        dev=[abs(a-b)/max(b,1e-9) for a,b,_,n in effs if b and b>0.05]
        stat['mc_vs_exact_maxrel']=float(max(dev)) if dev else None
    # ---- 3. tally chains
    nb=0; ucache={}
    for cd in dc:
        m=re.match(r'^F(\d+):P\s+(.*)$',cd,re.I)
        if not m or '<' not in m.group(2): continue
        body=m.group(2)
        # split top-level bins
        bins=[]; depth=0; cur=''
        for ch in body:
            if ch=='(':
                depth+=1
            cur+=ch
            if ch==')':
                depth-=1
                if depth==0: bins.append(cur.strip()); cur=''
        for b in bins:
            parts_=b[1:-1].split('<')
            lo,hi=parts_[0],parts_[-1]
            oc=[int(x) for x in re.findall(r'\d+',lo)]; cc=[int(x) for x in re.findall(r'\d+',hi)]
            lats={cells[c]['_lat'] for c in cc}
            if len(parts_)==3:
                ll={int(x) for x in re.findall(r'\d+',parts_[1])}
                if ll!=lats: msgs.append(f'tally F{m.group(1)} bin {b[:50]}: lattice list does not match containers')
            us={cells[c]['u'] for c in oc}
            for L in lats:
                if L not in ucache: ucache[L]=set(np.unique(latinfo[L][2]).tolist())
            allpresent=set().union(*[ucache[L] for L in lats])
            if not us<=allpresent: msgs.append(f'tally F{m.group(1)} bin {b[:60]}: organ not in lattices')
            nb+=1
    return dict(msgs=msgs,src=stat,tally_bins=nb,containers=len(conts))
if __name__=='__main__':
    for fn in sys.argv[1:]:
        r=check(fn)
        print(os.path.basename(fn),'OK' if not r['msgs'] else 'PROBLEM',r['src'],'bins',r['tally_bins'],'conts',r['containers'],r['msgs'][:5],flush=True)

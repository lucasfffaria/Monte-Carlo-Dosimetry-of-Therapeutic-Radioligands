"""Independent consistency checker for the MCNP6 decks."""
import re, sys, os, glob, numpy as np, json
def read_cards(fn, base):
    raw=open(fn,newline='').read()
    assert '\t' not in raw, 'TAB found'
    lines=raw.replace('\r\n','\n').split('\n')
    title=lines[0]; blocks=[[]]; i=1
    for l in lines[1:]:
        if l.strip()=='' :
            blocks.append([]); continue
        blocks[-1].append(l)
    blocks=[b for b in blocks]
    def cards(blk):
        out=[]
        for l in blk:
            if len(l)>80: raise ValueError(f'line>80: {l}')
            if re.match(r'^[cC]( |$)',l): continue
            ll=l.split('$')[0].rstrip()
            if ll.lower().startswith('read '):
                f=re.search(r'file=(\S+)',ll,re.I).group(1)
                sub=open(os.path.join(base,f),newline='').read().replace('\r\n','\n').split('\n')
                out+=[[x] for x in cards([x for x in sub if x.strip()!=''])]
                continue
            if ll.startswith('     ') and out:
                out[-1].append(ll.strip())
            else:
                out.append([ll.strip()])
        return [' '.join(x) for x in out]
    return title, cards(blocks[0]), cards(blocks[1]), cards(blocks[2])
_FILLCACHE={}
def expand_fill(text):
    if text in _FILLCACHE: return _FILLCACHE[text]
    toks=text.split()
    vals=[]; reps=[]
    for t in toks:
        if t[-1] in 'rR':
            reps[-1]+=int(t[:-1])
        else:
            vals.append(int(t)); reps.append(1)
    arr=np.repeat(np.array(vals,dtype=np.int64),np.array(reps))
    _FILLCACHE[text]=arr
    return arr
def parse(fn):
    base=os.path.dirname(fn)
    title,cc,sc,dc=read_cards(fn,base)
    cells={}; surfs={}; data=[]
    for c in cc:
        t=c.split()
        num=int(t[0]); mat=int(t[1]); k=2
        rho=None
        if mat!=0: rho=float(t[2]); k=3
        rest=' '.join(t[k:])
        params={}
        m=re.search(r'fill=([^ ]+)(\s*\(([^)]*)\))?',rest)
        geo=rest
        p_idx=min([rest.find(x) for x in ('u=','imp:','fill=','lat=') if rest.find(x)>=0])
        geo=rest[:p_idx].split(); par=rest[p_idx:]
        for kw in ('u','lat'):
            mm=re.search(kw+r'=(-?\d+)',par)
            if mm: params[kw]=int(mm.group(1))
        mm=re.search(r'imp:p=(\d+)',par); params['imp']=int(mm.group(1))
        if 'fill=' in par:
            fpart=par.split('fill=')[1]
            mm=re.match(r'(-?\d+):(-?\d+)\s+(-?\d+):(-?\d+)\s+(-?\d+):(-?\d+)\s+(.*)$',fpart)
            if mm:
                i1,i2,j1,j2,k1,k2=map(int,mm.groups()[:6])
                vals=expand_fill(mm.group(7))
                n=(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
                assert len(vals)==n, (num,'fill count',len(vals),n)
                params['fillarr']=(i1,i2,j1,j2,k1,k2,vals)
            else:
                mm=re.match(r'(\d+)\s*\(([^)]*)\)',fpart)
                params['fill']=int(mm.group(1)); tr=[float(x) for x in mm.group(2).split()]
                params['tr']=tr
        assert num not in cells, ('dup cell',num)
        cells[num]=dict(mat=mat,rho=rho,geo=geo,**params)
    for s in sc:
        t=s.split(); n=int(t[0]); assert n not in surfs, ('dup surf',n)
        surfs[n]=(t[1].lower(),[float(x) for x in t[2:]])
    return title,cells,surfs,dc
def surf_value(s,p):
    typ,a=s
    x,y,z=p[...,0],p[...,1],p[...,2]
    if typ=='px': return x-a[0]
    if typ=='py': return y-a[0]
    if typ=='pz': return z-a[0]
    if typ=='so': return np.sqrt(x*x+y*y+z*z)-a[0]
    if typ=='s': return np.sqrt((x-a[0])**2+(y-a[1])**2+(z-a[2])**2)-a[3]
    if typ=='rpp':
        inside=(x>a[0])&(x<a[1])&(y>a[2])&(y<a[3])&(z>a[4])&(z<a[5])
        return np.where(inside,-1.0,1.0)
    raise ValueError(typ)
_SV={}
def cell_contains(cell,surfs,p):
    """evaluate MCNP cell geometry (intersection=blank, union=':', parentheses) at points p"""
    expr=' '.join(cell['geo']).replace('(',' ( ').replace(')',' ) ').replace(':',' : ').split()
    pos=[0]
    def sv(tok):
        sn=int(tok)
        if abs(sn) not in _SV: _SV[abs(sn)]=surf_value(surfs[abs(sn)],p)<0
        ins=_SV[abs(sn)]
        return ins if sn<0 else ~ins
    def parse_union():
        r=parse_inter()
        while pos[0]<len(expr) and expr[pos[0]]==':':
            pos[0]+=1; r=r|parse_inter()
        return r
    def parse_inter():
        r=np.ones(p.shape[0],bool)
        while pos[0]<len(expr) and expr[pos[0]] not in (':',')'):
            t=expr[pos[0]]
            if t=='(':
                pos[0]+=1; q=parse_union(); assert expr[pos[0]]==')'; pos[0]+=1; r=r&q
            else:
                pos[0]+=1; r=r&sv(t)
        return r
    res=parse_union(); assert pos[0]==len(expr),('geometry parse',cell['geo'])
    return res
def check_deck(fn,npts=150000,verbose=True):
    _SV.clear()
    title,cells,surfs,dc=parse(fn)
    msgs=[]
    # materials
    mats=set(); dists={}; tallies={}
    for c in dc:
        t=c.split(); k=t[0].upper()
        if re.match(r'^M\d+$',k): mats.add(int(k[1:]))
        m=re.match(r'^(SI|SP|DS|SB)(\d+)$',k)
        if m:
            dists.setdefault(int(m.group(2)),{})[m.group(1)]=t[1:]
    used_mats={c['mat'] for c in cells.values() if c['mat']!=0}
    assert used_mats<=mats, ('missing materials',used_mats-mats)
    # surfaces used
    for n,c in cells.items():
        for tok in ' '.join(c['geo']).replace(':',' ').replace('(',' ').replace(')',' ').split():
            assert abs(int(tok)) in surfs, ('missing surface',n,tok)
    # universes
    univ={}
    for n,c in cells.items():
        if 'u' in c: univ.setdefault(c['u'],[]).append(n)
    for n,c in cells.items():
        if 'fill' in c: assert c['fill'] in univ, ('fill universe missing',n,c['fill'])
        if 'fillarr' in c:
            u=set(np.unique(c['fillarr'][6]).tolist())
            missing=[x for x in u if x not in univ]
            assert not missing, ('fill array universes missing',n,missing[:5])
    # level-0 cells overlap/gap check by random points in world sphere
    lvl0=[n for n,c in cells.items() if 'u' not in c]
    world=[s for s in surfs if surfs[s][0]=='so' and s!=99]
    R=max(surfs[s][1][0] for s in world)
    rng=np.random.default_rng(1)
    # sample: half uniform in world, half concentrated in the bounding region of objects
    pts=rng.uniform(-R,R,(npts,3)); pts=pts[np.linalg.norm(pts,axis=1)<R*0.999]
    boxes=[surfs[s][1] for s in surfs if surfs[s][0]=='rpp']
    if boxes:
        b=np.array(boxes); lo=b[:,[0,2,4]].min(0); hi=b[:,[1,3,5]].max(0)
        p2=rng.uniform(lo,hi,(npts,3)); pts=np.vstack([pts,p2])
    cnt=np.zeros(len(pts),int)
    for n in lvl0:
        cnt+=cell_contains(cells[n],surfs,pts)
    outside=np.linalg.norm(pts,axis=1)>=R
    bad=(cnt!=1)
    # points exactly on boundaries are measure-zero; any bad point is an error
    if bad.any():
        idx=np.flatnonzero(bad)[:5]
        msgs.append(f'GEOMETRY: {bad.sum()} of {len(pts)} points not in exactly one level-0 cell, e.g. {pts[idx]} counts {cnt[idx]}')
    # source checks
    sd=[c for c in dc if c.upper().startswith('SDEF')]
    nent=None
    if sd:
        s=sd[0].upper()
        if 'CEL=D1' in s:
            si1=dists[1]['SI']; assert si1[0].upper()=='L'
            ent=re.findall(r'\(([^)]*)\)',' '.join(si1[1:]))
            nent=len(ent)
            sp1=dists[1]['SP']; assert len(sp1)==nent,('SP1 count',len(sp1),nent)
            for d in (2,3,4,5):
                ds=dists[d]['DS']; assert ds[0].upper()=='S'; assert len(ds)-1==nent,('DS count',d,len(ds)-1,nent)
                for x in ds[1:]: assert int(x) in dists,('missing dist',x)
            # each path: organ cell in universe filling lattice in container
            for p in ent:
                parts=[x.strip() for x in p.split('<')]
                cont=int(parts[-1].lstrip('-')); lat=int(parts[-2])
                assert 'fill' in cells[cont] and cells[lat].get('u')==cells[cont]['fill'], ('path',p)
                if len(parts)==3:
                    oc=int(parts[0]); fa=cells[lat]['fillarr'][6]
                    assert cells[oc]['u'] in set(np.unique(fa).tolist()),('organ not in lattice',p)
    # tallies: SD counts
    fcards={}
    for c in dc:
        m=re.match(r'^F(\d+):P\s+(.*)$',c,re.I)
        if m: fcards[int(m.group(1))]=m.group(2)
    for n,body in fcards.items():
        bins=re.findall(r'\(\(?[^()]*\)?[^()]*\)',body)
        # robust: count top-level parenthesis groups
        depth=0; nb=0; plain=0; tok=''
        for ch in body:
            if ch=='(':
                if depth==0: nb+=1
                depth+=1
            elif ch==')': depth-=1
        if nb==0: nb=len(body.split())
        sdl=[c for c in dc if re.match(rf'^SD{n}\s',c,re.I)]
        if sdl:
            assert len(sdl[0].split())-1==nb, ('SD count',n,len(sdl[0].split())-1,nb)
        # check referenced cells exist
        for x in re.findall(r'\d+',body):
            assert int(x) in cells, ('tally cell missing',n,x)
        de=[c for c in dc if re.match(rf'^DE{n}\s',c,re.I)]; df=[c for c in dc if re.match(rf'^DF{n}\s',c,re.I)]
        if de: assert len(de[0].split())==len(df[0].split()),('DE/DF',n)
        fu=[c for c in dc if re.match(rf'^FU{n}\s',c,re.I)]
        if fu:
            for x in fu[0].split()[1:]: assert int(x) in dists,('FU dist missing',n,x)
    for d in dists:
        assert 1<=d<=999
        if 'SI' in dists[d] and 'SP' in dists[d]:
            si=dists[d]['SI']; sp=dists[d]['SP']
            opt=si[0].upper()
            body=' '.join(si[1:]) if opt in ('H','L','A','S') else ' '.join(si)
            nsi=len(re.findall(r'\([^)]*\)',body)) if '(' in body else len(body.split())
            nsp=len(sp)-1 if sp[0].upper() in ('D','C','V','W') else len(sp)
            assert nsi==nsp,('SI/SP count',d,nsi,nsp)
    return dict(title=title,ncells=len(cells),nsurf=len(surfs),nent=nent,ntally=len(fcards),msgs=msgs,npts=len(pts))
if __name__=='__main__':
    import paths
    fns=sys.argv[1:] or sorted(glob.glob(str(paths.INPUTS/'*.i')))
    for fn in fns:
        try:
            r=check_deck(fn)
            print(os.path.basename(fn),'OK' if not r['msgs'] else 'PROBLEM',{k:r[k] for k in ('ncells','nsurf','nent','ntally','npts')}, r['msgs'])
        except Exception as e:
            print(os.path.basename(fn),'ERROR',repr(e)[:300])

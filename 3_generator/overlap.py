"""Analytic pairwise overlap test between 'object' cells (inside of one RPP / union of RPPs / sphere)."""
import sys, glob, os, numpy as np, re
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import validate
def parse(fn):
    import tempfile, shutil
    txt=open(fn,newline='').read().replace('\r\n','\n').split('\n')
    txt=[l for l in txt if not l.lower().startswith('read ')]
    tmp=fn+'.noread.tmp'
    open(tmp,'w').write('\n'.join(txt))
    try:
        title,cc,sc,dc=validate.read_cards(tmp,os.path.dirname(fn))
    finally:
        os.remove(tmp)
    cells={}; surfs={}
    for c in cc:
        t=c.split(); num=int(t[0]); mat=int(t[1]); k=2 if mat==0 else 3
        rest=' '.join(t[k:])
        idx=[rest.find(x) for x in ('u=','imp:','fill=','lat=') if rest.find(x)>=0]
        p=min(idx); geo=rest[:p].split(); par=rest[p:]
        d=dict(geo=geo)
        if re.search(r'\bu=',par): d['u']=1
        cells[num]=d
    for s in sc:
        t=s.split(); surfs[int(t[0])]=(t[1].lower(),[float(x) for x in t[2:]])
    return title,cells,surfs,dc
def objects(cells,surfs):
    objs=[]
    for n,c in cells.items():
        if 'u' in c: continue
        g=' '.join(c['geo'])
        terms=[t.split() for t in g.split(':')]
        if all(len(t)==1 and t[0].startswith('-') for t in terms):
            shapes=[]
            for t in terms:
                s=surfs[abs(int(t[0]))]
                shapes.append(s)
            objs.append((n,shapes))
    return objs
def box_of(s):
    typ,a=s
    if typ=='rpp': return np.array(a)
    if typ=='s': return np.array([a[0]-a[3],a[0]+a[3],a[1]-a[3],a[1]+a[3],a[2]-a[3],a[2]+a[3]])
    if typ=='so': return np.array([-a[0],a[0],-a[0],a[0],-a[0],a[0]])
def overlap_vol(s1,s2):
    b1,b2=box_of(s1),box_of(s2)
    lo=np.maximum(b1[[0,2,4]],b2[[0,2,4]]); hi=np.minimum(b1[[1,3,5]],b2[[1,3,5]])
    d=hi-lo
    if (d<=1e-6).any(): return 0.0
    if s1[0]=='rpp' and s2[0]=='rpp': return float(np.prod(d))
    # sphere involved: sample
    rng=np.random.default_rng(0); p=rng.uniform(lo,hi,(20000,3))
    def inside(s,p):
        typ,a=s
        if typ=='rpp': return (p[:,0]>a[0])&(p[:,0]<a[1])&(p[:,1]>a[2])&(p[:,1]<a[3])&(p[:,2]>a[4])&(p[:,2]<a[5])
        if typ=='s': return np.linalg.norm(p-np.array(a[:3]),axis=1)<a[3]
    f=(inside(s1,p)&inside(s2,p)).mean()
    return f*np.prod(d)
def check(fn):
    title,cells,surfs,dc=parse(fn)
    objs=objects(cells,surfs)
    # exclude big 'region' boxes: cabin slabs etc are not single '-s' cells (they have outside tokens) -> already excluded
    bad=[]
    for i in range(len(objs)):
        for j in range(i+1,len(objs)):
            for s1 in objs[i][1]:
                for s2 in objs[j][1]:
                    v=overlap_vol(s1,s2)
                    if v>1e-6: bad.append((objs[i][0],objs[j][0],v))
    return len(objs),bad
if __name__=='__main__':
    for fn in sys.argv[1:]:
        n,bad=check(fn); print(os.path.basename(fn),'objects',n,'overlaps',bad[:10])

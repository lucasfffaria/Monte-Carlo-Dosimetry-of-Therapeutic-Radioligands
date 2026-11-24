"""Analytic pairwise overlap test between level-0 objects for v2 decks.
Objects: cells whose geometry is a union of '-surface' terms (slab containers, detectors, seats,
shell parts) and envelope cells (their union part; their holes are their own slab containers).
Pairs with the same owner (envelope n and its slab containers n*100+s) are skipped, but each slab
container must lie inside its envelope union."""
import sys, os, re, numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from overlap import parse, box_of, overlap_vol
def objects(cells,surfs):
    objs=[]
    for n,c in cells.items():
        if 'u' in c or n in (9998,9999): continue
        g=' '.join(c['geo'])
        m=re.match(r'^\(([^()]*)\)\s*(.*)$',g)
        if m: head,tail=m.group(1),m.group(2).split()
        else:
            toks=g.split(); head=toks[0] if toks else ''; tail=toks[1:]
            if ':' in g and not m: head,tail=g,[]
        terms=[t.strip() for t in head.split(':')]
        if not all(re.fullmatch(r'-\d+',t) for t in terms): continue
        if tail and not all(re.fullmatch(r'\d+',t) for t in tail): continue
        shapes=[surfs[int(t[1:])] for t in terms]
        owner=n if tail else (n//100 if n>=100000 else n)
        objs.append((n,shapes,owner,[int(t) for t in tail]))
    return objs
def check(fn):
    title,cells,surfs,dc=parse(fn)
    objs=objects(cells,surfs)
    bad=[]
    B=[np.array([box_of(s) for s in o[1]]) for o in objs]
    lo=np.array([b[:,[0,2,4]].min(0) for b in B]); hi=np.array([b[:,[1,3,5]].max(0) for b in B])
    def in_hole(i,j):
        # object i lies inside one of the holes of object j
        if not objs[j][3]: return False
        def inside(s,hb): return (box_of(s)[[0,2,4]]>=hb[[0,2,4]]-1e-9).all() and (box_of(s)[[1,3,5]]<=hb[[1,3,5]]+1e-9).all()
        hbs=[box_of(surfs[h]) for h in objs[j][3]]
        return all(any(inside(s,hb) for hb in hbs) for s in objs[i][1])
    for i in range(len(objs)):
        # candidate partners by bounding box
        cand=np.flatnonzero(((np.minimum(hi[i],hi)-np.maximum(lo[i],lo))>1e-6).all(1))
        for j in cand:
            if j<=i or objs[i][2]==objs[j][2]: continue
            if in_hole(i,j) or in_hole(j,i): continue
            for s1 in objs[i][1]:
                for s2 in objs[j][1]:
                    v=overlap_vol(s1,s2)
                    if v>1e-6: bad.append((objs[i][0],objs[j][0],round(v,6)))
    # slab containers inside their envelope union and listed as holes
    env={o[0]:o for o in objs if o[3]}
    miss=[]
    for o in objs:
        if o[0]>=100000:
            e=env.get(o[0]//100)
            if e is None: miss.append((o[0],'no envelope')); continue
            sb=box_of(o[1][0])
            ok=any((sb[[0,2,4]]>=box_of(s)[[0,2,4]]).all() and (sb[[1,3,5]]<=box_of(s)[[1,3,5]]).all() for s in e[1])
            if not ok: miss.append((o[0],'not inside envelope'))
            if int(cells[o[0]]['geo'][0][1:]) not in e[3]: miss.append((o[0],'not a hole of envelope'))
    return len(objs),bad,miss
if __name__=='__main__':
    for fn in sys.argv[1:]:
        n,bad,miss=check(fn)
        print(os.path.basename(fn),'objects',n,'overlaps',len(bad),bad[:6],'containment',miss[:6],flush=True)

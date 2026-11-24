import sys, os, numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from overlap import parse, box_of
for fn in sys.argv[1:]:
    title,cells,surfs,dc=parse(fn)
    bad=[]; nchk=0
    for n,c in cells.items():
        if 'u' in c: continue
        g=c['geo']
        if len(g)>1 and g[0].startswith('-') and ':' not in ' '.join(g):
            outer=surfs[abs(int(g[0]))]
            if outer[0] not in ('rpp','so'): continue
            ob=box_of(outer)
            for tok in g[1:]:
                s=surfs[abs(int(tok))]; b=box_of(s); nchk+=1
                if s[0]=='rpp' and outer[0]=='rpp':
                    inside=(b[0]>=ob[0]-1e-9)&(b[1]<=ob[1]+1e-9)&(b[2]>=ob[2]-1e-9)&(b[3]<=ob[3]+1e-9)&(b[4]>=ob[4]-1e-9)&(b[5]<=ob[5]+1e-9)
                    # objects listed in a region cell must lie entirely inside it (or entirely outside, e.g. world exclusions)
                    disjoint=(b[1]<=ob[0])|(b[0]>=ob[1])|(b[3]<=ob[2])|(b[2]>=ob[3])|(b[5]<=ob[4])|(b[4]>=ob[5])
                    if not (inside or disjoint): bad.append((n,tok,'partial'))
                elif outer[0]=='so':
                    r=outer[1][0]; corners=np.array([[b[i],b[2+j],b[4+k]] for i in (0,1) for j in (0,1) for k in (0,1)])
                    if (np.linalg.norm(corners,axis=1)>=r).any(): bad.append((n,tok,'outside world'))
                elif s[0]=='s' and outer[0]=='rpp':
                    if not ((b[0]>=ob[0])&(b[1]<=ob[1])&(b[2]>=ob[2])&(b[3]<=ob[3])&(b[4]>=ob[4])&(b[5]<=ob[5])): bad.append((n,tok,'sphere partial'))
    print(os.path.basename(fn),'checked',nchk,'problems',bad[:10])

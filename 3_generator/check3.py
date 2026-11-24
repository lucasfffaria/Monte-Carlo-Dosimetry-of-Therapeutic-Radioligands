"""v3 checks: energy-resolved source (D6 grid, TME=FCEL=D2 region tags, T0 bins, FT SCX 6) and tally sizes."""
import sys, os, re, numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import validate
def check(fn):
    title,cells,surfs,dc=validate.parse(fn)
    msgs=[]; info={}
    cards={}
    for c in dc:
        t=c.split(); cards.setdefault(t[0].upper(),[]).append(t[1:])
    sd=[c for c in dc if c.upper().startswith('SDEF')][0].upper()
    tallies={}
    for c in dc:
        m=re.match(r'^F(\d+):P\s+(.*)$',c,re.I)
        if m: tallies[int(m.group(1))]=m.group(2)
    def nbins(body):
        depth=0; nb=0
        for ch in body:
            if ch=='(':
                if depth==0: nb+=1
                depth+=1
            elif ch==')': depth-=1
        return nb if nb else len(body.split())
    if 'CEL=D1' in sd:
        if 'ERG=D6' not in sd or 'TME=FCEL=D2' not in sd: msgs.append('SDEF not energy-resolved')
        si6=cards['SI6'][0]; sp6=[float(x) for x in cards['SP6'][0]]
        E=[float(x) for x in si6[1:]]
        if si6[0].upper()!='L' or len(E)!=27 or len(sp6)!=27: msgs.append('D6 grid size')
        if abs(sum(sp6)-1)>1e-5: msgs.append(f'SP6 sum {sum(sp6)}')
        if min(E)<0.02 or max(E)<2.6: msgs.append('grid range')
        ds2=[int(x) for x in cards['DS2'][0][1:]]
        regs=sorted(set(ds2))
        if regs!=list(range(971,995)): msgs.append(f'DS2 region dists {regs[:3]}..{regs[-3:]} n={len(regs)}')
        for r in regs:
            si=cards[f'SI{r}'][0]
            if si[0].upper()!='L' or float(si[1])!=1000*(r-970): msgs.append(f'time tag {r} {si}')
        t0=[float(x) for x in cards['T0'][0]]
        if len(t0)!=24 or any(abs(t0[i]-1000*(i+1.5))>1e-6 for i in range(24)): msgs.append('T0 bins')
        for n in tallies:
            ft=cards.get(f'FT{n}')
            if not ft or ' '.join(ft[0]).upper()!='SCX 6': msgs.append(f'FT{n} not SCX 6: {ft}')
        if any(k.startswith('FU') for k in cards): msgs.append('FU card present')
        ntime=25; nuser=28
    else:
        ntime=1; nuser=len(cards['SI1'][0])-1+1
    # DE ranges
    for k,v in cards.items():
        if k.startswith('DE') and k[2:].isdigit():
            e=[float(x) for x in v[0][1:]]
            if max(e)<2.62: msgs.append(f'{k} max energy {max(e)}')
    total=0
    for n,body in tallies.items():
        nb=nbins(body)
        ne=len(cards[f'E{n}'][0])+1 if f'E{n}' in cards else 1
        total+=nb*ntime*nuser*ne
    info['tally_values']=total
    info['mctal_MB_est']=round(total*2*12/1e6,1)
    return msgs,info
if __name__=='__main__':
    for fn in sys.argv[1:]:
        msgs,info=check(fn)
        print(os.path.basename(fn),'OK' if not msgs else 'PROBLEM',info,msgs[:6],flush=True)

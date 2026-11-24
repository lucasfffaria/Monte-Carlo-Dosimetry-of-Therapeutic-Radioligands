"""Scene builders: verification (ICRP116), standing pairs (Article 1), bus (Article 2)."""
import numpy as np, os, json
from mcnpgen import *
import data as D
# ---------------------------------------------------------------- common geometry numbers
STD_TC={'AM':(26.42,15.27),'AF':(26.68,13.95)}      # trunk centre (local x,y) at abdomen level
def std_box(sex):
    ph=VARIANTS['TAMstd' if sex=='AM' else 'TAFstd'].ph
    return ph.extent()
def seated_info(sex):
    v=VARIANTS['TAMseat' if sex=='AM' else 'TAFseat']; ph=v.ph; m=ph.meta
    zs=m['z_seat_surface']; k=np.arange(ph.nk)+ph.k0; z=(k+0.5)*ph.dz
    below=(ph.arr>0)&(z[:,None,None]<zs); above=(ph.arr>0)&(z[:,None,None]>=zs)
    bb=ph.body_box(); bA=ph.body_box(above); bB=ph.body_box(below)
    return dict(zseat=zs,body=bb,A=bA,B=bB,xc=0.5*(bb[0]+bb[1]))
PATIENT_STD_XF=Xform((-STD_TC['AM'][0],-STD_TC['AM'][1],0.0))
def patient_standing(deck):
    P=Placement('patient','',None,7000) if False else None
    v=VARIANTS['PAMstd']
    P=Placement('PAT',v,PATIENT_STD_XF,7000)
    deck.place(P); return P
def target_std_placement(sex,where,gap,cell):
    """where in {'front','back','side'}; gap (cm) between bounding boxes along the facing axis"""
    v=VARIANTS['TAMstd' if sex=='AM' else 'TAFstd']
    e=v.ph.extent(); xc,yc=STD_TC[sex]
    pb=PATIENT_STD_XF.box(VARIANTS['PAMstd'].ph.extent())   # patient global box
    if where=='front':
        xf=Xform((xc, pb[2]-gap, 0.0),rot180=True)
    elif where=='back':
        xf=Xform((-xc, pb[3]+gap, 0.0))
    elif where=='side':
        xf=Xform((pb[1]+gap, -yc, 0.0))
    return Placement(f'{where.upper()}{sex}',v,xf,cell)
# ---------------------------------------------------------------- V: ICRP 116 verification
V_ENERGIES=[0.02,0.03,0.04,0.05,0.06,0.08,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.5,2.0,3.0]
def deck_V(sex,geom,nps,ctme,test=False):
    v=VARIANTS['TAMstd' if sex=='AM' else 'TAFstd']
    d=Deck(f'V {sex} {geom}: ICRP116 verification, broad parallel beam, ICRP110 {sex}',world_radius=300.0,world_void=True)
    d.notes=['Verification of the target dosimetry chain against ICRP 116 (E/phi, organ doses).',
             'Fluence = 1/A_beam per source photon; A_beam printed below. FT SCX 1 bins = photon energies (18).']
    T=Placement(f'{sex}',v,Xform((0,0,0)),7001)
    d.place(T)
    d.finish_world()
    e=v.ph.extent(); x0,x1,y0,y1,z0,z1=e
    if geom=='AP':   pos=('X=D3',f'Y={fpos(y0-0.5)}','Z=D4'); vec='0 1 0'; A=(x1-x0)*(z1-z0); rng=((x0,x1),(z0,z1))
    if geom=='PA':   pos=('X=D3',f'Y={fpos(y1+0.5)}','Z=D4'); vec='0 -1 0'; A=(x1-x0)*(z1-z0); rng=((x0,x1),(z0,z1))
    if geom=='LLAT': pos=(f'X={fpos(x1+0.5)}','Y=D3','Z=D4'); vec='-1 0 0'; A=(y1-y0)*(z1-z0); rng=((y0,y1),(z0,z1))
    if geom=='RLAT': pos=(f'X={fpos(x0-0.5)}','Y=D3','Z=D4'); vec='1 0 0'; A=(y1-y0)*(z1-z0); rng=((y0,y1),(z0,z1))
    d.notes.append(f'Beam area A = {A:.4f} cm2 (fluence per source photon = 1/A).')
    L=['c ============ SOURCE: broad parallel photon beam',
       f'SDEF PAR=P ERG=D1 {pos[0]} {pos[1]} {pos[2]} VEC={vec} DIR=1',
       'SI1 L '+' '.join(fsig(e,4) for e in V_ENERGIES),'SP1 '+' '.join('1' for e in V_ENERGIES),
       f'SI3 H {fpos(rng[0][0])} {fpos(rng[0][1])}','SP3 D 0 1',f'SI4 H {fpos(rng[1][0])} {fpos(rng[1][1])}','SP4 D 0 1']
    d.data+=L
    d.regions=[]
    target_tallies(d,[T],binning='SCX',scx_dist=1)
    run_cards(d,nps,ctme,test)
    d.meta=dict(kind='V',sex=sex,geom=geom,beam_area=A,energies=V_ENERGIES)
    return d
# ---------------------------------------------------------------- P0: patient alone + air detectors
def p0_detectors():
    pb=PATIENT_STD_XF.box(VARIANTS['PAMstd'].ph.extent())
    faces={'front':(0,-1,pb[2]),'back':(0,1,pb[3]),'left':(1,0,pb[1]),'right':(-1,0,pb[0])}
    dets=[]
    def add(lbl,c,r): dets.append((lbl,np.array(c,float),r))
    radius={10:3.0,30:4.0,50:5.0,100:8.0,200:8.0,300:8.0}
    for nm,(ux,uy,f) in faces.items():
        for g in (10,30,50,100,200,300):
            hs=(90,110,135) if g in (30,50,100,200) else (110,)
            for h in hs:
                if ux==0: c=(0.0,f+uy*g,h)
                else: c=(f+ux*g,0.0,h)
                add(f'{nm}_g{g}_z{h}',c,radius[g])
    for ang in (45,135,225,315):
        a=np.radians(ang)
        for R in (75,125,225):
            # angle measured from +x (left) counter-clockwise; front is -y (270 deg)
            add(f'diag{ang}_R{R}_z110',(R*np.cos(a),R*np.sin(a),110.0),6.0)
    for nm in ('front','left'):
        ux,uy,f=faces[nm]
        for h in (20,45,70,160):
            c=(0.0,f+uy*100,h) if ux==0 else (f+ux*100,0.0,h)
            add(f'{nm}_g100_z{h}',c,8.0)
    # overlap check
    for i in range(len(dets)):
        for j in range(i+1,len(dets)):
            dd=np.linalg.norm(dets[i][1]-dets[j][1])
            assert dd>dets[i][2]+dets[j][2]+0.5,('detector overlap',dets[i][0],dets[j][0])
    for lbl,c,r in dets:   # not inside patient box
        assert not (pb[0]-r<c[0]<pb[1]+r and pb[2]-r<c[1]<pb[3]+r and c[2]-r<pb[5]), lbl
    return dets
def deck_P0(nps,ctme,test=False):
    d=Deck('P0: 177Lu-PSMA patient (ICRP110 AM, standing) alone in air; H*(10)/Ka map',world_radius=600.0)
    d.notes=['Patient alone in air; air-sphere detectors at 10 cm - 3 m.',
             'Detectors: air spheres; positions/labels listed in comments and in the generator JSON.']
    P=patient_standing(d)
    dets=p0_detectors(); cells=[]; labels=[]
    for n,(lbl,c,r) in enumerate(dets):
        s=d.sph(c,r,lbl); cn=7101+n; d.det_vol[cn]=4/3*np.pi*r**3
        d.cell(cn,AIR_MAT,AIR_RHO,[f'-{s}'],comment=lbl,outside_tokens=[str(s)])
        cells.append(cn); labels.append(lbl)
    d.finish_world()
    patient_source(d,P,line_z=(2.0,174.0))
    spec=[cells[labels.index(x)] for x in ('front_g100_z110','left_g100_z110','back_g100_z110')]
    detector_tallies(d,cells,labels,spectrum_cells=spec)
    run_cards(d,nps,ctme,test)
    d.meta=dict(kind='P0',detectors=[(l,list(map(float,c)),r) for l,c,r in dets])
    return d
# ---------------------------------------------------------------- T: patient + 3 exposed persons
def deck_T(sex,gap_cm,nps,ctme,test=False):
    d=Deck(f'T {sex} gap {gap_cm:g} cm: patient AM standing + {sex} front/side/back',world_radius=600.0)
    d.notes=[f'Exposed {sex} phantoms facing the patient (FRONT, rotated 180 deg), beside the patient',
             f'(SIDE, left, same orientation) and behind (BACK, facing the patient back); gap = {gap_cm:g} cm',
             'between body bounding boxes along the facing axis.']
    P=patient_standing(d)
    targets=[target_std_placement(sex,w,gap_cm,7001+i) for i,w in enumerate(('front','side','back'))]
    for t in targets: d.place(t)
    d.finish_world()
    patient_source(d,P,line_z=(2.0,174.0))
    target_tallies(d,targets)
    run_cards(d,nps,ctme,test)
    d.meta=dict(kind='T',sex=sex,gap=gap_cm,targets=[(t.name,list(t.xf.o),t.xf.rot180,t.env) for t in targets])
    return d
# ---------------------------------------------------------------- bus
BUS=dict(module=50.5,aisle=35.0,pitch=80.0,seat_top=38.5,cushion_depth=42.0,cushion_w=45.0,cushion_t=10.0,
         pan_t=0.2,back_t=8.0,shell_t=1.0,back_h=65.0,half_in=120.5,roof=200.0,rows=range(0,7),occ_rows=range(1,6),patient_row=3,
         y_front=-460.5,y_rear=369.5)
def seat_x():
    m=BUS['module']; a=BUS['aisle']
    return {'LW':-(a/2+1.5*m),'LA':-(a/2+0.5*m),'RA':(a/2+0.5*m),'RW':(a/2+1.5*m)}
def yb(r): return (r-BUS['patient_row'])*BUS['pitch']
def seated_placement(sex,role,row,seat,cell):
    v=VARIANTS[('PAMseat' if role=='P' else ('TAMseat' if sex=='AM' else 'TAFseat'))]
    si=seated_info(sex)
    bb,bA,bB=si['body'],si['A'],si['B']
    # forward shift so that the shanks (box B) stay in front of the cushion front edge
    yback_g=min(yb(row), yb(row)-BUS['cushion_depth']+(bb[3]-bB[3])-0.05)
    o=(seat_x()[seat]-si['xc'], yback_g-bb[3], BUS['seat_top']-si['zseat'])
    xf=Xform(o)
    eps=1e-3
    boxesA=(bA[0]-eps,bA[1]+eps,bA[2]-eps,bA[3]+eps,bA[4]+eps,bA[5]+eps)
    boxesB=(bB[0]-eps,bB[1]+eps,bB[2]-eps,bB[3]+eps,bB[4]-eps,bB[5]-eps)
    P=Placement(f'{role}{sex}_r{row}{seat}',v,xf,cell,env_local=[boxesA,boxesB])
    P.seat=(row,seat)
    return P
def seat_parts(row,seat):
    """global RPP boxes and materials of one seat"""
    xs=seat_x()[seat]; w=BUS['cushion_w']/2; y=yb(row); st=BUS['seat_top']
    parts=[('cushion',(xs-w,xs+w,y-BUS['cushion_depth'],y,st-BUS['cushion_t'],st),206,0.05),
           ('pan',(xs-w,xs+w,y-BUS['cushion_depth'],y,st-BUS['cushion_t']-BUS['pan_t'],st-BUS['cushion_t']),201,7.85),
           ('back',(xs-w,xs+w,y,y+BUS['back_t'],st,st+BUS['back_h']),206,0.05),
           ('shell',(xs-w,xs+w,y+BUS['back_t'],y+BUS['back_t']+BUS['shell_t'],st,st+BUS['back_h']),207,0.90)]
    return parts
def bus_detectors(patient_seat):
    dets=[]
    for r in BUS['occ_rows']:
        for s,xs in seat_x().items():
            if (r,s)==(BUS['patient_row'],patient_seat): continue
            for h in (15,40,75):
                dets.append((f'r{r}{s}_h{h}',np.array([xs,yb(r)-13.5,BUS['seat_top']+h]),8.0))
    for r in (2,3,4):
        for z in (100,130,160):
            dets.append((f'aisle_r{r}_z{z}',np.array([0.0,yb(r)-40.0,float(z)]),8.0))
    return dets
def build_bus(d,objects_by_row,with_bus=True):
    """objects_by_row: row -> list of (surfaces tokens for 'outside', ) ; writes seats, slabs, shell"""
    hi=BUS['half_in']; R=BUS['roof']
    slabs=[]
    for r in BUS['rows']:
        y=yb(r); slab=(-hi,hi,y-BUS['pitch']+9.5,y+9.5,0.0,R)
        slabs.append((r,slab))
    y0=BUS['y_front']; y1=BUS['y_rear']
    cn=8001
    if with_bus:
        # seats
        for r in BUS['rows']:
            for s in seat_x():
                for nm,b,mat,rho in seat_parts(r,s):
                    sf=d.rpp(b,f'seat r{r}{s} {nm}')
                    d.cell(cn,mat,rho,[f'-{sf}'],comment=f'seat r{r}{s} {nm}',top=False); cn+=1
                    objects_by_row.setdefault(r,[]).append(str(sf))
    # cabin air slabs
    for r,slab in slabs:
        sf=d.rpp(slab,f'cabin slab row {r}')
        d.cell(cn,AIR_MAT,AIR_RHO,[f'-{sf}']+objects_by_row.get(r,[]),comment=f'cabin air row {r}',top=False); cn+=1
        d._slabs=getattr(d,'_slabs',[])+[sf]
    fr=d.rpp((-hi,hi,y0,slabs[0][1][2],0.0,R),'cabin front'); d.cell(cn,AIR_MAT,AIR_RHO,[f'-{fr}'],comment='cabin air front',top=False); cn+=1
    re_=d.rpp((-hi,hi,slabs[-1][1][3],y1,0.0,R),'cabin rear'); d.cell(cn,AIR_MAT,AIR_RHO,[f'-{re_}'],comment='cabin air rear',top=False); cn+=1
    if not with_bus:
        return cn
    # shell
    ho=hi+2.0
    shell=[('floor vinyl',(-ho,ho,y0-1,y1+1,-0.2,0.0),205,1.40),
           ('floor plywood',(-ho,ho,y0-1,y1+1,-2.0,-0.2),204,0.60),
           ('floor steel',(-ho,ho,y0-1,y1+1,-2.2,-2.0),201,7.85),
           ('roof panel',(-ho,ho,y0-1,y1+1,R,R+0.3),207,0.90),
           ('roof gap',(-ho,ho,y0-1,y1+1,R+0.3,R+4.3),AIR_MAT,AIR_RHO),
           ('roof skin',(-ho,ho,y0-1,y1+1,R+4.3,R+4.5),202,2.70),
           ('front wall',(-hi,hi,y0-1,y0,0.0,R),207,0.90),
           ('rear wall',(-hi,hi,y1,y1+1,0.0,R),207,0.90)]
    for sgn,side in ((-1,'L'),(1,'R')):
        def xr(a,b): return (min(sgn*a,sgn*b),max(sgn*a,sgn*b))
        shell+=[(f'{side} panel low',(*xr(hi,hi+0.5),y0-1,y1+1,0.0,80.0),207,0.90),
                (f'{side} window',(*xr(hi,hi+0.5),y0-1,y1+1,80.0,175.0),203,2.50),
                (f'{side} panel up',(*xr(hi,hi+0.5),y0-1,y1+1,175.0,R),207,0.90),
                (f'{side} gap',(*xr(hi+0.5,ho-0.2),y0-1,y1+1,0.0,R),AIR_MAT,AIR_RHO),
                (f'{side} skin',(*xr(ho-0.2,ho),y0-1,y1+1,0.0,R),202,2.70)]
    env=d.rpp((-ho,ho,y0-1,y1+1,-2.2,R+4.5),'bus envelope')
    shell_s=[]
    for nm,b,mat,rho in shell:
        sf=d.rpp(b,nm); shell_s.append(sf)
        d.cell(cn,mat,rho,[f'-{sf}'],comment=nm,top=False); cn+=1
    d.outside.append(str(env))
    # check that the shell + cabin pieces tile the envelope: envelope remainder (corners) -> air
    toks=[f'-{env}']+[str(s) for s in shell_s]+[str(s) for s in d._slabs]+[str(fr),str(re_)]
    d.cell(cn,AIR_MAT,AIR_RHO,toks,comment='envelope remainder',top=False); cn+=1
    return cn
def deck_bus(kind,patient_seat,sex,nps,ctme,test=False):
    """kind: 'D' detectors in bus, 'DAIR' detectors without bus, 'F' full occupancy (targets sex)"""
    tag={'RW':'window','RA':'aisle'}[patient_seat]
    if kind=='F':
        title=f'BUS-F {sex}: patient AM seated ({tag} seat), all other 19 seats occupied by {sex}'
    elif kind=='D':
        title=f'BUS-D: patient AM seated ({tag}) in empty bus; H*(10)/Ka at seats'
    else:
        title=f'BUS-DAIR: patient AM seated ({tag}) in air (no bus); same detectors'
    d=Deck(title,world_radius=700.0)
    d.notes=['Coach 2+2 (EMTU/AGER based): module 50.5 cm, aisle 35 cm, pitch 80 cm,',
             'seat height 38.5 cm, cushion 45x42 cm. Rows 1-5 occupied, patient row 3.']
    pr=BUS['patient_row']
    P=seated_placement('AM','P',pr,patient_seat,7000)
    objects={}
    d.place(P,top=(kind=='DAIR'))
    objects.setdefault(pr,[]).extend(str(s) for s in P.env_surfs)
    targets=[]; dets=[]
    if kind=='F':
        n=0
        for r in BUS['occ_rows']:
            for s in seat_x():
                if (r,s)==(pr,patient_seat): continue
                t=seated_placement(sex,'T',r,s,7001+n); n+=1
                d.place(t,top=False); targets.append(t)
                objects.setdefault(r,[]).extend(str(x) for x in t.env_surfs)
    else:
        dets=bus_detectors(patient_seat); cells=[]; labels=[]
        for n,(lbl,c,rr) in enumerate(dets):
            s=d.sph(c,rr,lbl); cn=7101+n; d.det_vol[cn]=4/3*np.pi*rr**3
            d.cell(cn,AIR_MAT,AIR_RHO,[f'-{s}'],comment=lbl,top=(kind=='DAIR'),outside_tokens=[str(s)])
            # which row slab?
            row=int(round(BUS['patient_row']+(c[1]+13.5)/BUS['pitch'])) if not lbl.startswith('aisle') else int(lbl.split('_r')[1][0])
            objects.setdefault(row,[]).append(str(s))
            cells.append(cn); labels.append(lbl)
    if kind!='DAIR':
        build_bus(d,objects,with_bus=True)
    d.finish_world()
    patient_source(d,P,line_z=(si_line:=(seated_info('AM')['zseat'],176.0)))
    if kind=='F': target_tallies(d,targets)
    else: detector_tallies(d,cells,labels,spectrum_cells=[cells[labels.index(x)] for x in labels if x.endswith('_h40') and ('r3' in x)][:3])
    run_cards(d,nps,ctme,test)
    d.meta=dict(kind='BUS-'+kind,patient_seat=patient_seat,sex=sex,
                targets=[(t.name,t.seat,list(t.xf.o),t.env) for t in targets],
                detectors=[(l,list(map(float,c)),r) for l,c,r in dets],
                patient=(P.name,list(P.xf.o),P.env))
    return d

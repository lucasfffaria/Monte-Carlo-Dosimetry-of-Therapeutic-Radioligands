import os, sys, json, time, shutil
import numpy as np
from mcnpgen import *
import scenes, paths
# output folder for the rebuilt decks (compare them with the official ones); can be given on the command line
OUT=sys.argv[1] if len(sys.argv)>1 else str(paths.DECKS_OUT)
RUNS=[ # name, builder, args, nps, ctme(min), phase  -- NPS set from the v4 test speeds (01/11/2026):
       # V 640-740, P0 385, T 275-285, B01-03 275-315, B04-07 110-126 M histories/h (20 threads)
 ('V01_AM_AP',  scenes.deck_V,('AM','AP'),  1e8,30,1),
 ('V02_AM_PA',  scenes.deck_V,('AM','PA'),  1e8,30,1),
 ('V03_AM_LLAT',scenes.deck_V,('AM','LLAT'),1e8,30,1),
 ('V04_AM_RLAT',scenes.deck_V,('AM','RLAT'),1e8,30,1),
 ('V05_AF_AP',  scenes.deck_V,('AF','AP'),  1e8,30,1),
 ('V06_AF_PA',  scenes.deck_V,('AF','PA'),  1e8,30,1),
 ('V07_AF_LLAT',scenes.deck_V,('AF','LLAT'),1e8,30,1),
 ('V08_AF_RLAT',scenes.deck_V,('AF','RLAT'),1e8,30,1),
 ('P0_paciente_em_pe',scenes.deck_P0,(),     1e9,300,2),
 ('T01_AM_g030',scenes.deck_T,('AM',30),     8e8,400,2),
 ('T02_AM_g050',scenes.deck_T,('AM',50),     8e8,400,2),
 ('T03_AM_g100',scenes.deck_T,('AM',100),    8e8,400,2),
 ('T04_AM_g200',scenes.deck_T,('AM',200),    8e8,400,2),
 ('T05_AF_g030',scenes.deck_T,('AF',30),     8e8,400,2),
 ('T06_AF_g050',scenes.deck_T,('AF',50),     8e8,400,2),
 ('T07_AF_g100',scenes.deck_T,('AF',100),    8e8,400,2),
 ('T08_AF_g200',scenes.deck_T,('AF',200),    8e8,400,2),
 ('B01_DAIR_janela',scenes.deck_bus,('DAIR','RW','AM'),8e8,400,3),
 ('B02_D_janela',   scenes.deck_bus,('D','RW','AM'),   8e8,400,3),
 ('B03_D_corredor', scenes.deck_bus,('D','RA','AM'),   8e8,400,3),
 ('B04_F_AM_janela',  scenes.deck_bus,('F','RW','AM'), 1e9,900,3),
 ('B05_F_AF_janela',  scenes.deck_bus,('F','RW','AF'), 1e9,1000,3),
 ('B06_F_AM_corredor',scenes.deck_bus,('F','RA','AM'), 1e9,900,3),
 ('B07_F_AF_corredor',scenes.deck_bus,('F','RA','AF'), 1e9,1000,3),
]
SPEED={'V':690e6,'P0':385e6,'TAM':284e6,'TAF':275e6,'BD':290e6,'BFAM':126e6,'BFAF':111e6}   # histories per hour (tests)
def expected_hours(name,nps):
    if name.startswith('V'): k='V'
    elif name.startswith('P0'): k='P0'
    elif name.startswith('T'): k='TAM' if '_AM_' in name else 'TAF'
    elif name[:3] in ('B01','B02','B03'): k='BD'
    else: k='BFAM' if '_AM_' in name else 'BFAF'
    return nps/SPEED[k]
TESTS=['V01_AM_AP','V05_AF_AP','P0_paciente_em_pe','T01_AM_g030','T05_AF_g030',
       'B01_DAIR_janela','B02_D_janela','B03_D_corredor','B04_F_AM_janela','B05_F_AF_janela']
DIAG=[('DX1_V01_sem_tallies','none'),('DX2_V01_formato_v3','omit'),('DX3_V01_por_fatia','perslab')]
TEST_ORDER=TESTS
def build(only=None):
    os.makedirs(OUT,exist_ok=True)
    meta={}
    for name,fn,args,nps,ctme,ph in RUNS:
        if only and name not in only: continue
        t=time.time()
        d=fn(*args,nps,ctme,False); d.write(f'{OUT}/{name}.i')
        if name in TESTS:
            dt=fn(*args,3e5 if name.startswith('V') else 1e6,20,True); dt.title=(name+'_TESTE '+dt.title)[:78]; dt.write(f'{OUT}/{name}_teste.i')
        m=getattr(d,'meta',{}); m['phase']=ph; m['nps']=nps; m['ctme']=ctme
        m['tally_index']=getattr(d,'tally_index',{})
        if hasattr(d,'src_entries'):
            ent=d.src_entries
            m['regions']=[(r[0],r[1]) for r in d.regions]
            tot=sum(e['p'] for e in ent)
            m['region_prob']={str(970+i+1):sum(e['p'] for e in ent if e['ed']==970+i+1)/tot for i in range(len(d.regions))}
            m['n_source_entries']=len(ent)
            m['source_stats']=d.src_stats
            m['energy_grid_MeV']=d.egrid; m['energy_prob']=d.pE
            m['region_time_tag_shakes']={r[0]:1000.0*(i+1) for i,r in enumerate(d.regions)}
        m['variants']=sorted(d.variants); m['expected_hours']=round(expected_hours(name,nps),2)
        m['placements']=[dict(name=P.name,variant=P.var.key,envelope_cell=P.envcell,slab_cells=P.slabcells,
                              origin=list(P.xf.o),rot180=P.xf.rot180,envelope=P.env) for P in d.placements]
        meta[name]=m
        print(f'{name:22s} {time.time()-t:6.1f}s  size {os.path.getsize(f"{OUT}/{name}.i")/1e3:8.1f} kB  vars {sorted(d.variants)}',flush=True)
    return meta
if __name__=='__main__':
    meta=build()
    def conv(o):
        if isinstance(o,(np.floating,)): return float(o)
        if isinstance(o,(np.integer,)): return int(o)
        if isinstance(o,np.ndarray): return o.tolist()
        return str(o)
    json.dump(meta,open(OUT+'/metadados_runs.json','w',newline='\n'),indent=1,default=conv)
    write_lattice_files(OUT)
    for k,v in VARIANTS.items():
        print('lattice',k,os.path.getsize(os.path.join(OUT,v.latfile()))/1e6,'MB')

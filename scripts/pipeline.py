"""Extended pipeline (v25). Reuses Marco's run_all.py functions; adds solar-control film runs,
night-only service level, demand-tariff (BTD) billing, PV-first counterfactual, evening load,
fixed critical load, and incremental economics. Output: results.json"""
import json, math, sys
from pathlib import Path
R=Path(__file__).resolve().parents[1]
src=(R/'scripts/run_all.py').read_text()
src=src[:src.index('all_city={}')].replace("ROOT = Path(__file__).resolve().parents[1]",f"ROOT = Path(r'{R}')").replace("cache=True","cache=False")
g={}; exec(src,g)
pd=g['pd']; np=g['np']; CITY=g['CITY']; A=g['ASSUMPTIONS']; CAPEX=g['CAPEX']; TAR=g['TARIFFS']
dispatch=g['dispatch']; read_load=g['read_load']; knee=g['knee_capacity']; outage=g['outage_test']; minres=g['minimum_reserve_for_48h']
FX=A['exchange_rd_per_usd']; r=A['discount_rate_real']; N=A['horizon_years']
DEM={'EDENORTE':(207.14,9.83,1310.10),'EDEESTE':(208.78,9.07,1280.88)}   # BTD, SIE-044-2026-TF, Apr-Jun 2026
EF=0.6277
def prof(code,case,svc):
    name={('E0','full'):f'{code}_E0',('E1','full'):f'{code}_E1_film',('E2','full'):f'{code}_E2_film',
          ('E0','night'):f'{code}_E0_night',('E1','night'):f'{code}_E1_night',('E2','night'):f'{code}_E2_night'}[(case,svc)]
    d=read_load(R/f'inputs/energyplus/hourly_profiles_v25/{name}_hourly.csv').rename(columns={'total_electric_kWh':'load_kWh'}); return d
mid=lambda a,b:(a+b)/2
K={'E1':mid(CAPEX['E1_low'],CAPEX['E1_high']),'E2':mid(CAPEX['E2_low'],CAPEX['E2_high'])}
PVK=CAPEX['pv_mid']; BATR=mid(CAPEX['battery_usd_kwh_low'],CAPEX['battery_usd_kwh_high'])
# ---------- billing ----------
def bts_bill(load_df,dist,imp=None,exp=None,gen_annual=0,ppa=0,accounts=1):
    t=TAR[dist]; mo=load_df.timestamp.dt.month.to_numpy(); L=load_df.load_kWh.to_numpy()
    imp=L if imp is None else imp; exp=np.zeros_like(L) if exp is None else exp
    tot=0; credit=0
    for m in range(1,13):
        sel=mo==m; ip=imp[sel].sum()/accounts; xp=exp[sel].sum()/accounts
        net=ip-xp-credit; billed=max(net,0); credit=max(-net,0)
        fixed=(t['fixed_high'] if L[sel].sum()/accounts>100 else t['fixed_low'])
        tot+=(g['energy_charge'](billed,t)+fixed)*accounts + exp[sel].sum()*A['network_charge_fraction']*t['r4']
    pay=min(credit,A['annual_credit_generation_limit']*gen_annual/accounts)*accounts*ppa/1000*FX
    return tot-pay
def btd_bill(load_df,dist,imp=None,exp=None,gen_annual=0,ppa=0):
    f,e,d=DEM[dist]; mo=load_df.timestamp.dt.month.to_numpy(); L=load_df.load_kWh.to_numpy()
    imp=L if imp is None else imp; exp=np.zeros_like(L) if exp is None else exp
    mmax=np.array([imp[mo==m].max() for m in range(1,13)]); top2=np.sort(mmax)[-2:].mean(); bd=np.maximum(mmax,top2)
    tot=0; credit=0
    for m in range(1,13):
        sel=mo==m; net=imp[sel].sum()-exp[sel].sum()-credit; billed=max(net,0); credit=max(-net,0)
        tot+=f+billed*e+bd[m-1]*d
    pay=min(credit,A['annual_credit_generation_limit']*gen_annual)*ppa/1000*FX
    return tot-pay, bd.mean()
def bills(code,loads,gen,cap,tariff,deg=1.0,floor=None):
    dist=CITY[code]['distributor']; ppa=CITY[code]['ppa_usd_mwh']; out={}
    for s in ('E0','E1','E2'):
        out[s]= bts_bill(loads[s],dist) if tariff=='flat' else (bts_bill(loads[s],dist,accounts=42) if tariff=='indiv' else btd_bill(loads[s],dist)[0])
    l2=loads['E2']; G=gen*deg
    for s,c in (('E3',0),('E4',cap)):
        rr=dispatch(G,l2.load_kWh.to_numpy(),c,**({'operating_floor_frac':floor} if (floor and c>0) else {}))
        if tariff=='flat': out[s]=bts_bill(l2,dist,rr['import'],rr['export'],G.sum(),ppa)
        elif tariff=='indiv': out[s]=bts_bill(l2,dist,rr['import'],rr['export'],G.sum(),ppa,accounts=42)
        else: out[s]=btd_bill(l2,dist,rr['import'],rr['export'],G.sum(),ppa)[0]
    return out
def npv_irr(cfs):
    npv=sum(c/(1+r)**i for i,c in enumerate(cfs))
    lo,hi=-0.9,5
    f=lambda x: sum(c/(1+x)**i for i,c in enumerate(cfs))
    if f(lo)*f(hi)>0: irr=float('nan')
    else:
        for _ in range(200):
            m=(lo+hi)/2
            if f(lo)*f(m)<=0: hi=m
            else: lo=m
        irr=(lo+hi)/2
    cum=0; pb=None
    for i,c in enumerate(cfs):
        cum+=c
        if cum>=0 and pb is None: pb=i
    return npv,100*irr,pb
def finance(code,loads,gen,cap,tariff,floor=None):
    # PV degradation applied year by year for E3/E4 (bill recomputed each 5 years and interpolated for speed)
    b0=bills(code,loads,gen,cap,tariff,1.0,floor); yrs=list(range(1,N+1))
    anchors={y:bills(code,loads,gen,cap,tariff,(1-A['pv_degradation'])**(y-1),floor) for y in (1,6,11,16,21,25)}
    def bill(s,y):
        if s in ('E0','E1','E2'): return b0[s]
        ks=sorted(anchors); lo=max(k for k in ks if k<=y); hi=min(k for k in ks if k>=y)
        if lo==hi: return anchors[lo][s]
        w=(y-lo)/(hi-lo); return anchors[lo][s]*(1-w)+anchors[hi][s]*w
    capex={'E1':K['E1'],'E2':K['E2'],'E3':K['E2']+PVK,'E4':K['E2']+PVK+cap*BATR}
    res={}
    for s in ('E1','E2','E3','E4'):
        cfs=[-capex[s]]; sav=[0]; om=[0]; rep=[0]
        for y in yrs:
            sv=(b0['E0']-bill(s,y))/FX; o_=(A['pv_om_fraction']*PVK if s in ('E3','E4') else 0); rp=(cap*BATR if (s=='E4' and y==A['battery_replacement_year']) else 0)
            cfs.append(sv-o_-rp); sav.append(sv); om.append(o_); rep.append(rp)
        n,i,p=npv_irr(cfs); res[s]={'K':capex[s],'sav1':sav[1],'NPV':n,'IRR':i,'PB':p,'sav':sav,'om':om,'rep':rep}
    res['bill1_RD']={s:b0[s] for s in b0}
    return res
out={}
for code,cfg in CITY.items():
    pv,_=g['pv_profile'](cfg['epw']); gen=pv.pv_ac_kWh.to_numpy(); o={'PV_MWh':gen.sum()/1000,'POA':pv.poa_global_Wh_m2.sum()/1000,'GHI':None}
    for svc in ('full','night'):
        loads={s:prof(code,s,svc) for s in ('E0','E1','E2')}
        q={}
        for s,d in loads.items():
            q[s]={'load_MWh':d.load_kWh.sum()/1000,'cool_MWh':d.cooling_electric_kWh.sum()/1000,'cool_th_MWh':d.cooling_thermal_kWh.sum()/1000,
                  'latent_MWh':d.latent_kWh.sum()/1000,'latent_pct':100*d.latent_kWh.sum()/d.cooling_thermal_kWh.sum(),
                  'peak_cool_kW':d.cooling_electric_kWh.max(),'peak_kW':d.load_kWh.max(),'kwh_apt_month':d.load_kWh.sum()/42/12,
                  'evening_kW':d.load_kWh[(d.timestamp.dt.hour>=18)&(d.timestamp.dt.hour<=22)].mean()}
        L2=loads['E2'].load_kWh.to_numpy(); L0=loads['E0'].load_kWh.to_numpy()
        cap,sweep=knee(gen,L2); r3=dispatch(gen,L2,0); r4=dispatch(gen,L2,cap)
        ev=((loads['E2'].timestamp.dt.hour>=18)&(loads['E2'].timestamp.dt.hour<=22)).to_numpy()
        cap0,_=knee(gen,L0); p3=dispatch(gen,L0,0); p4=dispatch(gen,L0,cap0); p4b=dispatch(gen,L0,cap)
        q['E3']={'SCR':100*r3['self'].sum()/gen.sum(),'SSR':100*r3['self'].sum()/L2.sum(),'imp_MWh':r3['import'].sum()/1000,'exp_MWh':r3['export'].sum()/1000,'evening_imp_kW':r3['import'][ev].mean()}
        q['E4']={'cap':cap,'SCR':100*r4['self'].sum()/gen.sum(),'SSR':100*r4['self'].sum()/L2.sum(),'imp_MWh':r4['import'].sum()/1000,'exp_MWh':r4['export'].sum()/1000,'loss_MWh':r4['losses'].sum()/1000,'evening_imp_kW':r4['import'][ev].mean(),
                 'sweep':sweep[['battery_kWh','SSR_pct']].round(2).values.tolist()}
        q['PVfirst']={'SSR_E0PV':100*p3['self'].sum()/L0.sum(),'imp_E0PV':p3['import'].sum()/1000,'exp_E0PV':p3['export'].sum()/1000,'cap0':cap0,'SSR_E0PVbat':100*p4['self'].sum()/L0.sum(),'imp_E0PVbat':p4['import'].sum()/1000,'SSR_E0PVsamebat':100*p4b['self'].sum()/L0.sum(),'imp_E0PVsamebat':p4b['import'].sum()/1000}
        if svc=='full':
            crit=A['critical_load_fraction']*L2.mean(); n48=outage(gen,cap,r4['soc'],crit,48); best=minres(gen,L2,cap,crit)
            critfix=6.0
            nf=outage(gen,cap,r4['soc'],critfix,48); bestf=minres(gen,L2,cap,critfix)
            q['E5']={'crit_kW':crit,'aut_h':cap*0.9*math.sqrt(0.9)/crit,'ord_pct':n48['success_pct'],'mean_served':n48['mean_served_pct'],'reserve':best[0]*100,
                     'fixed6kW_ord_pct':nf['success_pct'],'fixed6kW_reserve':(bestf[0]*100 if bestf else None),
                     'dur':{dd:outage(gen,cap,r4['soc'],crit,dd)['success_pct'] for dd in (4,12,24,72,168)}}
            # cost of the reserve: E4 dispatch with reserve -> SSR and imports
            rr=dispatch(gen,L2,cap,operating_floor_frac=best[0]); q['E5']['SSR_reserve']=100*rr['self'].sum()/L2.sum(); q['E5']['imp_reserve_MWh']=rr['import'].sum()/1000
        q['fin']={t:finance(code,loads,gen,cap,t) for t in ('BTD','flat','indiv')}
        if svc=='full': q['fin_reserve']=finance(code,loads,gen,cap,'BTD',floor=best[0])
        q['btd_billed_kW']={s:btd_bill(loads[s],CITY[code]['distributor'])[1] for s in ('E0','E1','E2')}
        q['hourly_E0']=loads['E0'].load_kWh.round(3).tolist() if False else None
        o[svc]=q
        print(code,svc,'done',flush=True)
    out[code]=o
json.dump(out,open(R/'results/results.json','w'),indent=1,default=float)

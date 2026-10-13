import json, numpy as np
from pathlib import Path
R=str(Path(__file__).resolve().parents[1]/'results')+'/'
r=json.load(open(R+'results.json'))
SITES=['MC','SD','ST']; NAME={'MC':'Montecristi','SD':'Santo Domingo','ST':'Santiago'}
A=2943.0; EF=0.6277; FX=60.5; rr=0.08; N=25; a=(1-(1+rr)**-N)/rr; NA=42
DIST={'MC':'EDENORTE','SD':'EDEESTE','ST':'EDENORTE'}
TOP={'MC':14.04,'SD':13.26,'ST':14.04}; BTDE={'MC':9.83,'SD':9.07,'ST':9.83}
VENT={'MC':((229091,228642),(133417,131131)),'SD':((216801,211210),(131344,119844)),'ST':((164625,140875),(98877,75167))}
CAPR={'Passive':(40084,71069),'Active':(115728-40084,207197-71069),'PV':(95950,131300)}
def npv(c,r_=rr): return sum(x/(1+r_)**i for i,x in enumerate(c))
def irr(c):
    lo,hi=-0.9,5
    if npv(c,lo)*npv(c,hi)>0: return float('nan')
    for _ in range(200):
        m=(lo+hi)/2
        if npv(c,lo)*npv(c,m)<=0: hi=m
        else: lo=m
    return (lo+hi)/2
out={}
for s in SITES:
    o=r[s]; F=o['full']; Nn=o['night']; fin=F['fin']['BTD']; d={}
    E={k:F[k] for k in ('E0','E1','E2')}
    d['cool']={k:E[k]['cool_MWh'] for k in E}; d['eui']={k:E[k]['cool_MWh']*1000/A for k in E}
    d['red']={'E1':1-E['E1']['cool_MWh']/E['E0']['cool_MWh'],'E2':1-E['E2']['cool_MWh']/E['E0']['cool_MWh']}
    d['E1toE2']=1-E['E2']['cool_MWh']/E['E1']['cool_MWh']; d['resid']=1-E['E2']['cool_MWh']/(E['E1']['cool_MWh']*0.7)
    d['envequip']=1-E['E1']['cool_MWh']*0.7/E['E0']['cool_MWh']
    d['lat']={k:E[k]['latent_pct'] for k in E}; d['latMWh']={k:E[k]['latent_MWh'] for k in E}
    d['peak']={k:E[k]['peak_cool_kW'] for k in E}; d['peakred']=1-E['E2']['peak_cool_kW']/E['E0']['peak_cool_kW']
    d['load']={k:E[k]['load_MWh'] for k in E}; d['nc']={k:E[k]['load_MWh']-E[k]['cool_MWh'] for k in E}
    d['kwh_apt']={k:E[k]['kwh_apt_month'] for k in E}
    d['PV']=o['PV_MWh']; d['yield']=o['PV_MWh']*1000/101; d['POA']=o['POA']
    d['E3']=F['E3']; d['E4']={k:v for k,v in F['E4'].items() if k!='sweep'}; d['E5']=F['E5']
    d['ratio0']=o['PV_MWh']/E['E0']['load_MWh']; d['ratio2']=o['PV_MWh']/E['E2']['load_MWh']
    imp=[E['E0']['load_MWh'],E['E1']['load_MWh'],E['E2']['load_MWh'],F['E3']['imp_MWh'],F['E4']['imp_MWh']]
    d['imports']=imp; d['impchange']=imp[4]/imp[0]-1; d['effshare']=(imp[0]-imp[2])/(imp[0]-imp[4]); d['impcut_apt']=(imp[0]-imp[4])/NA
    d['evening']={k:E[k]['evening_kW'] for k in E}; d['evening']['E3']=F['E3']['evening_imp_kW']; d['evening']['E4']=F['E4']['evening_imp_kW']
    d['PVfirst']=F['PVfirst']
    d['co2_pv']=o['PV_MWh']*EF; d['co2_eff']=(E['E0']['load_MWh']-E['E2']['load_MWh'])*EF; d['co2_tot']=d['co2_pv']+d['co2_eff']
    d['unz0']=(o['PV_MWh']-d['nc']['E0'])/E['E0']['cool_MWh']; d['unz2']=(o['PV_MWh']-d['nc']['E2'])/E['E2']['cool_MWh']
    d['vent']={'E1':VENT[s][0][0]/VENT[s][0][1]-1,'E2':VENT[s][1][0]/VENT[s][1][1]-1}
    # night service
    d['night']={'cool':{k:Nn[k]['cool_MWh'] for k in ('E0','E1','E2')},'load':{k:Nn[k]['load_MWh'] for k in ('E0','E1','E2')},
                'kwh_apt':{k:Nn[k]['kwh_apt_month'] for k in ('E0','E1','E2')},'peak':{k:Nn[k]['peak_cool_kW'] for k in ('E0','E1','E2')},
                'lat':{k:Nn[k]['latent_pct'] for k in ('E0','E1','E2')},
                'ratio0':o['PV_MWh']/Nn['E0']['load_MWh'],'ratio2':o['PV_MWh']/Nn['E2']['load_MWh'],
                'E3':Nn['E3'],'E4':{k:v for k,v in Nn['E4'].items() if k!='sweep'},'PVfirst':Nn['PVfirst'],
                'imports':[Nn['E0']['load_MWh'],Nn['E1']['load_MWh'],Nn['E2']['load_MWh'],Nn['E3']['imp_MWh'],Nn['E4']['imp_MWh']],
                'evening':{'E0':Nn['E0']['evening_kW'],'E2':Nn['E2']['evening_kW'],'E4':Nn['E4']['evening_imp_kW']}}
    # finance: cumulative & incremental
    def cfs(f,sc): return [ -f[sc]['K'] ] + [f[sc]['sav'][y]-f[sc]['om'][y]-f[sc]['rep'][y] for y in range(1,N+1)]
    def block(f):
        cum={sc:{'K':f[sc]['K'],'sav':f[sc]['sav1'],'NPV':f[sc]['NPV'],'IRR':f[sc]['IRR'],'PB':f[sc]['PB']} for sc in ('E1','E2','E3','E4')}
        prev=None; inc={}
        for sc,lab in (('E1','Passive'),('E2','Active'),('E3','PV'),('E4','Battery')):
            c=cfs(f,sc); p=[0]*(N+1) if prev is None else cfs(f,prev)
            ic=[x-y for x,y in zip(c,p)]
            dK=f[sc]['K']-(f[prev]['K'] if prev else 0); dS=f[sc]['sav1']-(f[prev]['sav1'] if prev else 0)
            om=f[sc]['om'][1]-(f[prev]['om'][1] if prev else 0)
            inc[lab]={'dK':dK,'dS':dS,'dNPV':npv(ic),'IRR':100*irr(ic),'PBs':dK/(dS-om) if dS-om>0 else float('nan'),'cf':ic}
            prev=sc
        return cum,inc
    d['fin']={}
    for t in ('BTD','flat','indiv'):
        cum,inc=block(F['fin'][t]); d['fin'][t]={'cum':cum,'inc':{k:{kk:vv for kk,vv in v.items() if kk!='cf'} for k,v in inc.items()}}
        if t=='BTD': incB=inc
    cumN,incN=block(Nn['fin']['BTD']); d['fin']['night_BTD']={'cum':cumN,'inc':{k:{kk:vv for kk,vv in v.items() if kk!='cf'} for k,v in incN.items()}}
    d['bill1']=F['fin']['BTD']['bill1_RD']; d['billed_kW']=F['btd_billed_kW']
    # resilience cost
    dn=incB['Battery']['dNPV']; d['res_cost']=-dn/(12*NA*a)
    nres=F['fin_reserve']['E4']['NPV']-F['fin']['BTD']['E3']['NPV']; d['res_cost_reserve']=-nres/(12*NA*a)
    # CCE / LCOE
    dE1=(E['E0']['load_MWh']-E['E1']['load_MWh'])*1000; dE2=(E['E1']['load_MWh']-E['E2']['load_MWh'])*1000
    d['cce']={'Passive':incB['Passive']['dK']/(a*dE1),'Active':incB['Active']['dK']/(a*dE2)}; d['dE']={'Passive':dE1/1000,'Active':dE2/1000}
    pvk=(95950+131300)/2; t_=np.arange(1,N+1)
    d['lcoe']=(pvk+sum(0.01*pvk/(1+rr)**t_))/sum(o['PV_MWh']*1000*0.995**(t_-1)/(1+rr)**t_)
    # switching values on incremental BTD cash flows
    sw={}
    for lab,v in incB.items():
        c=v['cf']; f=F['fin']['BTD']; sc={'Passive':'E1','Active':'E2','PV':'E3','Battery':'E4'}[lab]; pv_={'Passive':None,'Active':'E1','PV':'E2','Battery':'E3'}[lab]
        S_=[f[sc]['sav'][y]-(f[pv_]['sav'][y] if pv_ else 0) for y in range(N+1)]
        rest=[c[y]-S_[y] for y in range(N+1)]   # capex, om, replacement (negative)
        PVS=npv(S_); PVR=npv(rest)
        xs=-PVR/PVS if PVS>0 else float('nan')     # savings multiplier at NPV=0
        K0=-c[0]; other=npv(rest[1:]+[0]) if False else PVR+K0
        km=(PVS+other)/K0 if K0>0 else float('nan')
        sw[lab]={'IRR':100*irr(c),'rate_mult':xs,'cost_mult':km,'rate_usd':xs*BTDE[s]/FX}
        # tornado
        lo,hi=CAPR.get(lab,(None,None))
        if lab=='Battery': cap=F['E4']['cap']; lo,hi=250*cap,400*cap
        mid=(lo+hi)/2
        def npv_mod(kmult=1.0,r_=rr,g=0.0):
            cc=[c[0]*kmult]+[S_[y]*(1+g)**(y-1)+rest[y]*(kmult if (lab=='Battery' and rest[y]<-1000) else 1) for y in range(1,N+1)]
            return npv(cc,r_)
        sw[lab]['tornado']={'inv_hi':npv_mod(hi/mid),'inv_lo':npv_mod(lo/mid),'dr10':npv_mod(r_=0.10),'dr6':npv_mod(r_=0.06),'gm2':npv_mod(g=-0.02),'gp2':npv_mod(g=0.02),
                            'flat':d['fin']['flat']['inc'][lab]['dNPV'],'indiv':d['fin']['indiv']['inc'][lab]['dNPV'],'night':d['fin']['night_BTD']['inc'][lab]['dNPV'],'central':v['dNPV']}
    d['switch']=sw
    out[s]=d
json.dump(out,open(R+'derived.json','w'),indent=1,default=float)
def P(x,n=1): return round(x,n)
for s in SITES:
    d=out[s]; print('=====',s)
    print('cool',{k:P(v) for k,v in d['cool'].items()},'eui',{k:P(v) for k,v in d['eui'].items()},'red',{k:P(100*v) for k,v in d['red'].items()},'E1toE2',P(100*d['E1toE2']),'resid',P(100*d['resid']),'envequip',P(100*d['envequip']))
    print('lat',{k:P(v) for k,v in d['lat'].items()},'peak',{k:P(v) for k,v in d['peak'].items()},'peakred',P(100*d['peakred']),'load',{k:P(v) for k,v in d['load'].items()},'nc',{k:P(v) for k,v in d['nc'].items()},'kwh/apt',{k:P(v,0) for k,v in d['kwh_apt'].items()})
    print('PV',P(d['PV']),'yield',P(d['yield'],0),'ratio0/2',P(d['ratio0'],2),P(d['ratio2'],2),'E3',{k:P(v) for k,v in d['E3'].items()},'E4',{k:P(v) for k,v in d['E4'].items()})
    print('E5',d['E5'])
    print('imports',[P(x) for x in d['imports']],'chg',P(100*d['impchange']),'effshare',P(100*d['effshare']),'cut/apt',P(d['impcut_apt'],2),'evening',{k:P(v) for k,v in d['evening'].items()})
    print('PVfirst',{k:P(v) for k,v in d['PVfirst'].items()},'co2',P(d['co2_pv']),P(d['co2_eff']),P(d['co2_tot']),'unz',P(100*d['unz0']),P(100*d['unz2']),'vent',{k:P(100*v) for k,v in d['vent'].items()})
    print('night',json.dumps({k:(P(v,2) if isinstance(v,float) else v) for k,v in d['night'].items() if k not in ('E3','E4','PVfirst')},default=lambda x:round(x,1)))
    print('   nightE3',{k:P(v) for k,v in d['night']['E3'].items()},'E4',{k:P(v) for k,v in d['night']['E4'].items()},'PVfirst',{k:P(v) for k,v in d['night']['PVfirst'].items()})
    for t in ('BTD','flat','indiv','night_BTD'):
        print('  ',t,'cum',{k:(round(v['K']),round(v['sav']),round(v['NPV']),P(v['IRR']),v['PB']) for k,v in d['fin'][t]['cum'].items()})
        print('  ',t,'inc',{k:(round(v['dK']),round(v['dS']),round(v['dNPV']),P(v['IRR']),P(v['PBs'])) for k,v in d['fin'][t]['inc'].items()})
    print('bill1 kRD',{k:round(v/1000) for k,v in d['bill1'].items()},'billed kW',{k:P(v) for k,v in d['billed_kW'].items()})
    print('res cost',P(d['res_cost']),P(d['res_cost_reserve']),'cce',{k:P(v,3) for k,v in d['cce'].items()},'dE',{k:P(v) for k,v in d['dE'].items()},'lcoe',P(d['lcoe'],3))
    for k,v in d['switch'].items(): print('  sw',k,{kk:(P(vv,3) if not isinstance(vv,dict) else {a_:round(b_/1000,1) for a_,b_ in vv.items()}) for kk,vv in v.items()})

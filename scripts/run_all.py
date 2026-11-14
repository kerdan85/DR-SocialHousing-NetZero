from __future__ import annotations
from pathlib import Path
import json, math, hashlib, shutil, sys
import numpy as np
import pandas as pd
import pvlib
import matplotlib.pyplot as plt
from numba import njit

ROOT = Path(__file__).resolve().parents[1]
EP_BASE = ROOT / 'inputs/energyplus/validation_package'
OUT = ROOT / 'outputs'
FIG = ROOT / 'figures'
for p in [OUT/'hourly',OUT/'monthly',OUT/'summaries',FIG]: p.mkdir(parents=True,exist_ok=True)

CITY = {
    'MC': {
        'name': 'Las Matas - Montecristi', 'distributor': 'EDENORTE', 'zone': 'Cibao Noroeste',
        'epw': EP_BASE/'climate/DOM_MC_Monte.Cristi.784510_TMYx.epw',
        'profiles': {s: EP_BASE/f'hourly_profiles/MC_{s}_final046_hourly.csv' for s in ['E0','E1','E2']},
        'ppa_usd_mwh': 95.92,
    },
    'SD': {
        'name': 'Las Matas - Santo Domingo Este', 'distributor': 'EDEESTE', 'zone': 'Ozama o Metropolitana',
        'epw': EP_BASE/'climate/DOM_SD_Las.Americas.Intl.AP.784850_TMYx.epw',
        'profiles': {s: EP_BASE/f'hourly_profiles/SD_{s}_final046_hourly.csv' for s in ['E0','E1','E2']},
        'ppa_usd_mwh': 117.51,
    },
    'ST': {
        'name': 'Las Matas - Santiago', 'distributor': 'EDENORTE', 'zone': 'Cibao Norte',
        'epw': EP_BASE/'climate/DOM_ST_Santiago-Cibao.Intl.AP.784600_TMYx.epw',
        'profiles': {s: EP_BASE/f'hourly_profiles/ST_{s}_final046_hourly.csv' for s in ['E0','E1','E2']},
        'ppa_usd_mwh': 99.56,
    },
}

TARIFFS = {
    'EDENORTE': {'fixed_low':40.33,'fixed_high':126.81,'r1':5.97,'r2':8.51,'r3':13.83,'r4':14.04},
    'EDEESTE':  {'fixed_low':41.34,'fixed_high':127.83,'r1':6.17,'r2':8.71,'r3':13.04,'r4':13.26},
    'EDESUR':   {'fixed_low':42.10,'fixed_high':128.59,'r1':6.05,'r2':8.59,'r3':12.89,'r4':13.09},
}

ASSUMPTIONS = {
    'pv_kwp':101.0, 'roof_area_m2':505.0, 'module_efficiency':0.20, 'noct_c':45.0,
    'bos_factor':0.86, 'battery_roundtrip_eff':0.90, 'battery_min_soc':0.10,
    'battery_c_rate':0.50, 'critical_load_fraction':0.30, 'outage_hours':48,
    'exchange_rd_per_usd':60.50, 'discount_rate_real':0.08, 'horizon_years':25,
    'pv_degradation':0.005, 'pv_om_fraction':0.01, 'battery_replacement_year':13,
    'network_charge_fraction':0.25, 'annual_credit_generation_limit':0.50,
    'apartments':42,
}
CAPEX = {
    'roof_low':26168.0,'roof_high':44009.0,
    'solar_control_low':13916.0,'solar_control_high':27060.0,
    'hvac_low':40937.0,'hvac_high':57291.0,
    'led_low':17361.0,'led_high':41667.0,
    'dehum_low':17346.0,'dehum_high':37170.0,
    'pv_low':95950.0,'pv_high':131300.0,
    'battery_usd_kwh_low':250.0,'battery_usd_kwh_high':400.0,
}

def midpoint(a,b): return (a+b)/2
CAPEX['E1_low']=CAPEX['roof_low']+CAPEX['solar_control_low']
CAPEX['E1_high']=CAPEX['roof_high']+CAPEX['solar_control_high']
CAPEX['E2_low']=CAPEX['E1_low']+CAPEX['hvac_low']+CAPEX['led_low']+CAPEX['dehum_low']
CAPEX['E2_high']=CAPEX['E1_high']+CAPEX['hvac_high']+CAPEX['led_high']+CAPEX['dehum_high']
CAPEX['pv_mid']=midpoint(CAPEX['pv_low'],CAPEX['pv_high'])

def read_load(path: Path) -> pd.DataFrame:
    df=pd.read_csv(path)
    dt=pd.to_datetime(dict(year=df.year,month=df.month,day=df.day,hour=(df.hour_ending-1)%24))
    df.insert(0,'timestamp',dt)
    return df

def pv_profile(epw: Path) -> tuple[pd.DataFrame,dict]:
    data,meta=pvlib.iotools.read_epw(epw,coerce_year=2021)
    sp=pvlib.solarposition.get_solarposition(data.index,meta['latitude'],meta['longitude'],altitude=meta['altitude'],temperature=data['temp_air'])
    poa=pvlib.irradiance.get_total_irradiance(
        surface_tilt=0,surface_azimuth=180,solar_zenith=sp.apparent_zenith,solar_azimuth=sp.azimuth,
        dni=data.dni,ghi=data.ghi,dhi=data.dhi,albedo=0.2,model='haydavies',dni_extra=pvlib.irradiance.get_extra_radiation(data.index))
    mods=pvlib.pvsystem.retrieve_sam('CECMod')
    module_name='Canadian_Solar_Inc__CS1U_400MS'
    mod=mods[module_name]
    tc=pvlib.temperature.noct_sam(poa.poa_global,data.temp_air,data.wind_speed,noct=ASSUMPTIONS['noct_c'],module_efficiency=ASSUMPTIONS['module_efficiency'],transmittance_absorptance=0.9,array_height=1,mount_standoff=4)
    IL,I0,Rs,Rsh,nV=pvlib.pvsystem.calcparams_desoto(poa.poa_global,tc,alpha_sc=mod.alpha_sc,a_ref=mod.a_ref,I_L_ref=mod.I_L_ref,I_o_ref=mod.I_o_ref,R_sh_ref=mod.R_sh_ref,R_s=mod.R_s)
    sd=pvlib.pvsystem.singlediode(IL,I0,Rs,Rsh,nV,method='lambertw')
    ac=(sd.p_mp.fillna(0).clip(lower=0)*(ASSUMPTIONS['pv_kwp']*1000/mod.STC)*ASSUMPTIONS['bos_factor'])/1000
    out=pd.DataFrame({'timestamp':data.index.tz_localize(None),'pv_ac_kWh':ac.to_numpy(),'poa_global_Wh_m2':poa.poa_global.to_numpy(),'temp_air_C':data.temp_air.to_numpy(),'wind_speed_m_s':data.wind_speed.to_numpy()})
    details={'module':module_name,'module_stc_W':float(mod.STC),'latitude':meta['latitude'],'longitude':meta['longitude'],'altitude_m':meta['altitude']}
    return out,details

@njit(cache=True)
def _dispatch_core(gen,load,cap_kwh,operating_floor_frac,physical_min_frac,repeat,eta_rt,c_rate):
    n=len(gen)
    direct=np.zeros(n*repeat); bat=np.zeros(n*repeat); charge=np.zeros(n*repeat); socarr=np.zeros(n*repeat)
    imp=np.zeros(n*repeat); exp=np.zeros(n*repeat); selfe=np.zeros(n*repeat); losses=np.zeros(n*repeat)
    eta_c=math.sqrt(eta_rt); eta_d=eta_c; pmax=cap_kwh*c_rate
    floor=max(operating_floor_frac,physical_min_frac)*cap_kwh
    soc=0.5*cap_kwh
    for i in range(n*repeat):
        g=gen[i % n]; l=load[i % n]
        di=g if g<l else l
        surplus=g-l if g>l else 0.0; deficit=l-g if l>g else 0.0
        charge_in=surplus
        if charge_in>pmax: charge_in=pmax
        room=(cap_kwh-soc)/eta_c
        if charge_in>room: charge_in=room
        if charge_in<0: charge_in=0.0
        soc+=charge_in*eta_c
        discharge_out=deficit
        if discharge_out>pmax: discharge_out=pmax
        avail=(soc-floor)*eta_d
        if discharge_out>avail: discharge_out=avail
        if discharge_out<0: discharge_out=0.0
        soc-=discharge_out/eta_d
        direct[i]=di;bat[i]=discharge_out;charge[i]=charge_in;socarr[i]=soc
        imp[i]=deficit-discharge_out;exp[i]=surplus-charge_in;selfe[i]=di+discharge_out
        losses[i]=charge_in*(1-eta_c)+(discharge_out/eta_d-discharge_out)
    k=(repeat-1)*n
    return direct[k:],bat[k:],charge[k:],socarr[k:],imp[k:],exp[k:],selfe[k:],losses[k:]

def dispatch(gen,load,cap_kwh,operating_floor_frac=0.10,physical_min_frac=0.10,repeat=4,return_hourly=True):
    gen=np.asarray(gen,float);load=np.asarray(load,float);n=len(gen)
    if cap_kwh<=0:
        direct=np.minimum(gen,load); imp=np.maximum(load-gen,0); exp=np.maximum(gen-load,0)
        return {'direct':direct,'battery_to_load':np.zeros(n),'charge':np.zeros(n),'soc':np.zeros(n),'import':imp,'export':exp,'self':direct,'losses':np.zeros(n)}
    vals=_dispatch_core(gen,load,float(cap_kwh),float(operating_floor_frac),float(physical_min_frac),int(repeat),float(ASSUMPTIONS['battery_roundtrip_eff']),float(ASSUMPTIONS['battery_c_rate']))
    keys=['direct','battery_to_load','charge','soc','import','export','self','losses']
    return dict(zip(keys,vals))

def knee_capacity(gen,load):
    caps=np.arange(0,401,10,dtype=float);ssr=[];scr=[]
    for cap in caps:
        r=dispatch(gen,load,cap)
        self_e=r['self'].sum();ssr.append(100*self_e/load.sum());scr.append(100*self_e/gen.sum())
    x=(caps-caps.min())/(caps.max()-caps.min());y=(np.array(ssr)-min(ssr))/(max(ssr)-min(ssr));idx=int(np.argmax(y-x))
    sweep=pd.DataFrame({'battery_kWh':caps,'SSR_pct':ssr,'SCR_pct':scr,'distance_to_chord':y-x})
    return float(caps[idx]),sweep

@njit(cache=True)
def _outage_core(gen,cap_kwh,init_socs,critical_kw,duration,eta_rt,c_rate,min_soc):
    n=len(gen);eta_c=math.sqrt(eta_rt);eta_d=eta_c;pmax=cap_kwh*c_rate;physmin=cap_kwh*min_soc
    unmet=np.zeros(n)
    for start in range(n):
        soc=init_socs[start-1] if start>0 else init_socs[n-1]
        u=0.0
        for j in range(duration):
            g=gen[(start+j)%n];l=critical_kw
            surplus=g-l if g>l else 0.0;deficit=l-g if l>g else 0.0
            ch=surplus
            if ch>pmax: ch=pmax
            room=(cap_kwh-soc)/eta_c
            if ch>room: ch=room
            if ch<0: ch=0.0
            soc+=ch*eta_c
            dis=deficit
            if dis>pmax: dis=pmax
            avail=(soc-physmin)*eta_d
            if dis>avail: dis=avail
            if dis<0: dis=0.0
            soc-=dis/eta_d
            u+=deficit-dis
        unmet[start]=u
    return unmet

def outage_test(gen,cap_kwh,init_socs,critical_kw,duration=48):
    gen=np.asarray(gen,float);init_socs=np.asarray(init_socs,float)
    unmet=_outage_core(gen,float(cap_kwh),init_socs,float(critical_kw),int(duration),float(ASSUMPTIONS['battery_roundtrip_eff']),float(ASSUMPTIONS['battery_c_rate']),float(ASSUMPTIONS['battery_min_soc']))
    total=critical_kw*duration
    return {'success_pct':100*np.mean(unmet<1e-6),'mean_served_pct':100*(1-unmet.mean()/total),'minimum_served_pct':100*(1-unmet.max()/total),'mean_ENS_kWh':unmet.mean(),'max_ENS_kWh':unmet.max()}

def minimum_reserve_for_48h(gen,load,cap,critical):
    lo=int(ASSUMPTIONS['battery_min_soc']*100);best=None
    for pct in range(lo,101):
        r=dispatch(gen,load,cap,operating_floor_frac=pct/100)
        m=outage_test(gen,cap,r['soc'],critical,48)
        if m['success_pct']>=100-1e-9:
            best=(pct/100,r,m);break
    return best

def energy_charge(kwh,t):
    kwh=max(float(kwh),0)
    if kwh<=200: return kwh*t['r1']
    if kwh<=300: return 200*t['r1']+(kwh-200)*t['r2']
    if kwh<701: return 200*t['r1']+100*t['r2']+(kwh-300)*t['r3']
    return kwh*t['r4']

def monthly_bill_no_pv(load_df,distributor,accounts=1):
    t=TARIFFS[distributor]; records=[]
    for month,g in load_df.groupby(load_df.timestamp.dt.month):
        total=g.load_kWh.sum(); per=total/accounts
        fixed=(t['fixed_high'] if per>100 else t['fixed_low'])*accounts
        bill=energy_charge(per,t)*accounts+fixed
        records.append({'month':month,'gross_load_kWh':total,'grid_import_kWh':total,'grid_export_kWh':0,'net_billed_kWh':total,'credit_end_kWh':0,'energy_charge_RD':bill-fixed,'fixed_charge_RD':fixed,'network_charge_RD':0,'annual_credit_payment_RD':0.0,'total_bill_RD':bill})
    return pd.DataFrame(records)

def monthly_bill_net_meter(hourly,distributor,ppa_usd_mwh,accounts=1):
    t=TARIFFS[distributor];credit_per=0.0; records=[]
    annual_gen=hourly.pv_ac_kWh.sum()/accounts
    for month,g in hourly.groupby(hourly.timestamp.dt.month):
        gross_load=g.load_kWh.sum();gross_imp=g.grid_import_kWh.sum();gross_exp=g.grid_export_kWh.sum()
        imp_per=gross_imp/accounts;exp_per=gross_exp/accounts
        net=imp_per-exp_per-credit_per
        billed=max(net,0);credit_per=max(-net,0)
        fixed=(t['fixed_high'] if (g.load_kWh.sum()/accounts)>100 else t['fixed_low'])*accounts
        echarge=energy_charge(billed,t)*accounts
        # BTS1 is monomial. For the master meter, gross monthly withdrawal is >701 kWh; r4 is used as transparent proxy for the assigned energy charge.
        network=gross_exp*ASSUMPTIONS['network_charge_fraction']*t['r4']
        records.append({'month':month,'gross_load_kWh':gross_load,'grid_import_kWh':gross_imp,'grid_export_kWh':gross_exp,'net_billed_kWh':billed*accounts,'credit_end_kWh':credit_per*accounts,'energy_charge_RD':echarge,'fixed_charge_RD':fixed,'network_charge_RD':network,'annual_credit_payment_RD':0.0,'total_bill_RD':echarge+fixed+network})
    df=pd.DataFrame(records)
    allowed_credit_per=min(credit_per,ASSUMPTIONS['annual_credit_generation_limit']*annual_gen)
    payment=allowed_credit_per*accounts*(ppa_usd_mwh/1000)*ASSUMPTIONS['exchange_rd_per_usd']
    df.loc[df.index[-1],'annual_credit_payment_RD']=payment
    df.loc[df.index[-1],'total_bill_RD']-=payment
    return df

def npv(rate,cfs): return sum(cf/((1+rate)**i) for i,cf in enumerate(cfs))
def irr(cfs):
    lo,hi=-0.99,10.0
    flo=npv(lo,cfs);fhi=npv(hi,cfs)
    if flo*fhi>0:return float('nan')
    for _ in range(300):
        mid=(lo+hi)/2;fm=npv(mid,cfs)
        if abs(fm)<1e-8:return mid
        if flo*fm<=0:hi=mid;fhi=fm
        else:lo=mid;flo=fm
    return (lo+hi)/2

def payback(cfs,discount=0):
    cum=0
    for y,cf in enumerate(cfs):
        cum+=cf/((1+discount)**y)
        if cum>=0:return y
    return None

def finance(city_code,scenario,base_bill,loads,pv_base,cap,capex_level='mid',tariff_escalation=0,discount=None,tax_credit=False):
    discount=ASSUMPTIONS['discount_rate_real'] if discount is None else discount
    if capex_level=='low':
        e1=CAPEX['E1_low'];e2=CAPEX['E2_low'];pv_cap=CAPEX['pv_low'];bat_rate=CAPEX['battery_usd_kwh_low']
    elif capex_level=='high':
        e1=CAPEX['E1_high'];e2=CAPEX['E2_high'];pv_cap=CAPEX['pv_high'];bat_rate=CAPEX['battery_usd_kwh_high']
    else:
        e1=midpoint(CAPEX['E1_low'],CAPEX['E1_high']);e2=midpoint(CAPEX['E2_low'],CAPEX['E2_high']);pv_cap=CAPEX['pv_mid'];bat_rate=midpoint(CAPEX['battery_usd_kwh_low'],CAPEX['battery_usd_kwh_high'])
    capex={'E1':e1,'E2':e2,'E3':e2+pv_cap,'E4':e2+pv_cap+cap*bat_rate}[scenario]
    eligible=(pv_cap+(cap*bat_rate if scenario=='E4' else 0)) if scenario in ['E3','E4'] else 0
    if tax_credit:capex-=0.40*eligible
    cfs=[-capex]; year_bills=[]
    for y in range(1,ASSUMPTIONS['horizon_years']+1):
        escal=(1+tariff_escalation)**(y-1)
        if scenario in ['E1','E2']:
            bill=loads[scenario]['bill_RD']*escal
        else:
            deg=(1-ASSUMPTIONS['pv_degradation'])**(y-1)
            h=loads['E2']['hourly'].copy();g=pv_base*deg;l=h.load_kWh.to_numpy()
            r=dispatch(g,l,0 if scenario=='E3' else cap)
            h['pv_ac_kWh']=g;h['grid_import_kWh']=r['import'];h['grid_export_kWh']=r['export']
            b=monthly_bill_net_meter(h,CITY[city_code]['distributor'],CITY[city_code]['ppa_usd_mwh'],accounts=1)
            bill=b.total_bill_RD.sum()*escal
        year_bills.append(bill)
        savings=(base_bill*escal-bill)/ASSUMPTIONS['exchange_rd_per_usd']
        om=(ASSUMPTIONS['pv_om_fraction']*pv_cap if scenario in ['E3','E4'] else 0)
        cf=savings-om
        if scenario=='E4' and y==ASSUMPTIONS['battery_replacement_year']:cf-=cap*bat_rate
        cfs.append(cf)
    return {'capex_USD':capex,'cashflows_USD':cfs,'year1_savings_USD':cfs[1]+(ASSUMPTIONS['pv_om_fraction']*pv_cap if scenario in ['E3','E4'] else 0),'NPV_USD':npv(discount,cfs),'IRR_pct':100*irr(cfs),'PB_years':payback(cfs,0),'PBd_years':payback(cfs,discount),'year_bills_RD':year_bills}

def lcoe(production,pv_capex):
    pv_cost=pv_capex+sum(ASSUMPTIONS['pv_om_fraction']*pv_capex/(1+ASSUMPTIONS['discount_rate_real'])**y for y in range(1,ASSUMPTIONS['horizon_years']+1))
    pv_e=sum(production*((1-ASSUMPTIONS['pv_degradation'])**(y-1))/(1+ASSUMPTIONS['discount_rate_real'])**y for y in range(1,ASSUMPTIONS['horizon_years']+1))
    return pv_cost/pv_e

all_city={};e3_rows=[];e4_rows=[];e5_rows=[];fin_rows=[];bill_rows=[];monthly_all=[];sensitivity=[]
module_details={}
for code,cfg in CITY.items():
    loads={}
    for s,path in cfg['profiles'].items():
        d=read_load(path);d=d.rename(columns={'total_electric_kWh':'load_kWh'});loads[s]={'hourly':d,'annual_kWh':d.load_kWh.sum()}
        b=monthly_bill_no_pv(d,cfg['distributor'],1);loads[s]['billing']=b;loads[s]['bill_RD']=b.total_bill_RD.sum()
        b2=monthly_bill_no_pv(d,cfg['distributor'],ASSUMPTIONS['apartments']);loads[s]['bill_individual_RD']=b2.total_bill_RD.sum()
        b.assign(city=code,scenario=s,architecture='master_meter').to_csv(OUT/f'monthly/{code}_{s}_billing_master.csv',index=False)
    pv,md=pv_profile(cfg['epw']);module_details[code]=md
    e2=loads['E2']['hourly'].copy();gen=pv.pv_ac_kWh.to_numpy();load=e2.load_kWh.to_numpy();prod=gen.sum();dem=load.sum();direct=np.minimum(gen,load)
    h3=e2[['timestamp','load_kWh']].copy();h3['pv_ac_kWh']=gen;h3['direct_pv_kWh']=direct;h3['battery_to_load_kWh']=0.;h3['battery_charge_kWh']=0.;h3['soc_kWh']=0.;h3['grid_import_kWh']=np.maximum(load-gen,0);h3['grid_export_kWh']=np.maximum(gen-load,0);h3.to_csv(OUT/f'hourly/{code}_E3_hourly.csv',index=False)
    bill3=monthly_bill_net_meter(h3,cfg['distributor'],cfg['ppa_usd_mwh'],1);bill3.assign(city=code,scenario='E3',architecture='master_meter').to_csv(OUT/f'monthly/{code}_E3_billing_master.csv',index=False)
    bill3i=monthly_bill_net_meter(h3,cfg['distributor'],cfg['ppa_usd_mwh'],ASSUMPTIONS['apartments'])
    cap,sweep=knee_capacity(gen,load);sweep.to_csv(OUT/f'summaries/{code}_battery_sweep.csv',index=False)
    r4=dispatch(gen,load,cap);h4=e2[['timestamp','load_kWh']].copy();h4['pv_ac_kWh']=gen;h4['direct_pv_kWh']=r4['direct'];h4['battery_to_load_kWh']=r4['battery_to_load'];h4['battery_charge_kWh']=r4['charge'];h4['soc_kWh']=r4['soc'];h4['grid_import_kWh']=r4['import'];h4['grid_export_kWh']=r4['export'];h4['battery_losses_kWh']=r4['losses'];h4.to_csv(OUT/f'hourly/{code}_E4_hourly.csv',index=False)
    bill4=monthly_bill_net_meter(h4,cfg['distributor'],cfg['ppa_usd_mwh'],1);bill4.assign(city=code,scenario='E4',architecture='master_meter').to_csv(OUT/f'monthly/{code}_E4_billing_master.csv',index=False)
    bill4i=monthly_bill_net_meter(h4,cfg['distributor'],cfg['ppa_usd_mwh'],ASSUMPTIONS['apartments'])
    poa_kwh=pv.poa_global_Wh_m2.sum()/1000; pr=(prod/ASSUMPTIONS['pv_kwp'])/poa_kwh
    e3_rows.append({'city':code,'name':cfg['name'],'capacity_kWp':ASSUMPTIONS['pv_kwp'],'production_MWh':prod/1000,'specific_yield_kWh_kWp':prod/ASSUMPTIONS['pv_kwp'],'POA_kWh_m2':poa_kwh,'PR':pr,'E2_total_load_MWh':dem/1000,'SCR_pct':100*direct.sum()/prod,'SSR_pct':100*direct.sum()/dem,'annual_PV_load_ratio_pct':100*prod/dem,'grid_import_MWh':h3.grid_import_kWh.sum()/1000,'grid_export_MWh':h3.grid_export_kWh.sum()/1000,'LCOE_USD_kWh':lcoe(prod,CAPEX['pv_mid'])})
    e4_rows.append({'city':code,'name':cfg['name'],'battery_kWh':cap,'kWh_per_apartment':cap/ASSUMPTIONS['apartments'],'SSR_without_pct':100*direct.sum()/dem,'SSR_with_pct':100*r4['self'].sum()/dem,'SCR_with_pct':100*r4['self'].sum()/prod,'grid_import_MWh':r4['import'].sum()/1000,'grid_export_MWh':r4['export'].sum()/1000,'battery_losses_MWh':r4['losses'].sum()/1000})
    critical=ASSUMPTIONS['critical_load_fraction']*load.mean();normal48=outage_test(gen,cap,r4['soc'],critical,48);best=minimum_reserve_for_48h(gen,load,cap,critical);reserve,bestres,best48=best
    duration_rows=[]
    for dur in [12,24,48,72]:
        n=outage_test(gen,cap,r4['soc'],critical,dur);q=outage_test(gen,cap,bestres['soc'],critical,dur)
        duration_rows.append({'city':code,'duration_h':dur,'normal_success_pct':n['success_pct'],'resilient_success_pct':q['success_pct']})
    pd.DataFrame(duration_rows).to_csv(OUT/f'summaries/{code}_outage_duration_curve.csv',index=False)
    autonomy=cap*(1-ASSUMPTIONS['battery_min_soc'])*math.sqrt(ASSUMPTIONS['battery_roundtrip_eff'])/critical
    e5_rows.append({'city':code,'name':cfg['name'],'critical_load_kW':critical,'battery_kWh':cap,'autonomy_without_sun_h':autonomy,'normal_48h_success_pct':normal48['success_pct'],'normal_mean_served_pct':normal48['mean_served_pct'],'normal_minimum_served_pct':normal48['minimum_served_pct'],'normal_mean_ENS_kWh':normal48['mean_ENS_kWh'],'resilience_reserve_floor_pct':reserve*100,'resilient_48h_success_pct':best48['success_pct']})
    # billing and finance
    scenario_bills={'E0':loads['E0']['bill_RD'],'E1':loads['E1']['bill_RD'],'E2':loads['E2']['bill_RD'],'E3':bill3.total_bill_RD.sum(),'E4':bill4.total_bill_RD.sum()}
    scenario_bills_i={'E0':loads['E0']['bill_individual_RD'],'E1':loads['E1']['bill_individual_RD'],'E2':loads['E2']['bill_individual_RD'],'E3':bill3i.total_bill_RD.sum(),'E4':bill4i.total_bill_RD.sum()}
    for s in ['E0','E1','E2','E3','E4']:
        bill_rows.append({'city':code,'scenario':s,'distributor':cfg['distributor'],'annual_bill_master_RD':scenario_bills[s],'annual_bill_42_equal_accounts_RD':scenario_bills_i[s],'master_savings_vs_E0_RD':scenario_bills['E0']-scenario_bills[s],'individual_savings_vs_E0_RD':scenario_bills_i['E0']-scenario_bills_i[s]})
    loads['E3']={'bill_RD':scenario_bills['E3']};loads['E4']={'bill_RD':scenario_bills['E4']}
    for s in ['E1','E2','E3','E4']:
        f=finance(code,s,scenario_bills['E0'],loads,gen,cap)
        fin_rows.append({'city':code,'name':cfg['name'],'scenario':s,'investment_USD':f['capex_USD'],'year1_gross_savings_USD':f['year1_savings_USD'],'NPV_USD':f['NPV_USD'],'IRR_pct':f['IRR_pct'],'PB_years':f['PB_years'],'PBd_years':f['PBd_years']})
    # sensitivity E4: capex/discount/tariff/tax incentive
    for label,level,dr,esc,tax in [('Adverse','high',0.10,-0.02,False),('Central','mid',0.08,0.00,False),('Favorable','low',0.06,0.02,True)]:
        f=finance(code,'E4',scenario_bills['E0'],loads,gen,cap,capex_level=level,tariff_escalation=esc,discount=dr,tax_credit=tax)
        sensitivity.append({'city':code,'case':label,'capex_level':level,'discount_rate_pct':dr*100,'real_tariff_escalation_pct':esc*100,'tax_credit_40pct_eligible':tax,'NPV_USD':f['NPV_USD'],'IRR_pct':f['IRR_pct'],'investment_USD':f['capex_USD']})
    all_city[code]={'capacity_kWh':cap,'scenario_bills_RD':scenario_bills,'module_details':md}

# save summaries
for name,rows in [('E3_summary',e3_rows),('E4_summary',e4_rows),('E5_summary',e5_rows),('financial_summary',fin_rows),('billing_summary',bill_rows),('E4_sensitivity',sensitivity)]:
    pd.DataFrame(rows).to_csv(OUT/f'summaries/{name}.csv',index=False)
with open(OUT/'summaries/model_results.json','w',encoding='utf-8') as f:json.dump({'assumptions':ASSUMPTIONS,'capex':CAPEX,'tariffs':TARIFFS,'cities':all_city,'E3':e3_rows,'E4':e4_rows,'E5':e5_rows,'finance':fin_rows,'billing':bill_rows,'sensitivity':sensitivity,'module_details':module_details},f,ensure_ascii=False,indent=2)

# figures (no explicit colors to keep portable/default)
e3=pd.DataFrame(e3_rows);e4=pd.DataFrame(e4_rows);e5=pd.DataFrame(e5_rows);fin=pd.DataFrame(fin_rows)
labels=['Montecristi','Santo Domingo','Santiago']
fig,ax=plt.subplots(figsize=(8,5));x=np.arange(3);w=.35;ax.bar(x-w/2,e3.SCR_pct,w,label='SCR');ax.bar(x+w/2,e3.SSR_pct,w,label='SSR');ax.set_xticks(x,labels);ax.set_ylabel('%');ax.set_title('E3: autoconsumo y autosuficiencia sin batería');ax.legend();fig.tight_layout();fig.savefig(FIG/'Figure_E3_SCR_SSR.png',dpi=300);fig.savefig(FIG/'Figure_E3_SCR_SSR.svg');plt.close(fig)
fig,ax=plt.subplots(figsize=(8,5));ax.bar(x-w/2,e4.SSR_without_pct,w,label='Sin batería');ax.bar(x+w/2,e4.SSR_with_pct,w,label='Con batería');ax.set_xticks(x,labels);ax.set_ylabel('SSR (%)');ax.set_title('E4: efecto del almacenamiento');ax.legend();fig.tight_layout();fig.savefig(FIG/'Figure_E4_SSR.png',dpi=300);fig.savefig(FIG/'Figure_E4_SSR.svg');plt.close(fig)
fig,ax=plt.subplots(figsize=(8,5));ax.bar(x-w/2,e5.normal_48h_success_pct,w,label='Autoconsumo');ax.bar(x+w/2,e5.resilient_48h_success_pct,w,label='Reserva resiliente');ax.set_xticks(x,labels);ax.set_ylabel('Cortes de 48 h superados (%)');ax.set_ylim(0,105);ax.set_title('E5: resiliencia ante cortes de 48 horas');ax.legend();fig.tight_layout();fig.savefig(FIG/'Figure_E5_48h.png',dpi=300);fig.savefig(FIG/'Figure_E5_48h.svg');plt.close(fig)
fig,ax=plt.subplots(figsize=(9,5));pivot=fin.pivot(index='name',columns='scenario',values='NPV_USD').loc[[CITY[c]['name'] for c in ['MC','SD','ST']]];pivot.index=labels;pivot.plot(kind='bar',ax=ax);ax.set_ylabel('VAN (US$)');ax.set_title('Valor actual neto por escenario y provincia');ax.legend(title='Escenario');fig.tight_layout();fig.savefig(FIG/'Figure_financial_NPV.png',dpi=300);fig.savefig(FIG/'Figure_financial_NPV.svg');plt.close(fig)

# checksums
files=[p for p in ROOT.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt']
with open(ROOT/'SHA256SUMS.txt','w') as f:
    for p in sorted(files):
        h=hashlib.sha256(p.read_bytes()).hexdigest();f.write(f'{h}  {p.relative_to(ROOT)}\n')
print(pd.DataFrame(e3_rows).to_string(index=False))
print(pd.DataFrame(e4_rows).to_string(index=False))
print(pd.DataFrame(e5_rows).to_string(index=False))
print(pd.DataFrame(fin_rows).to_string(index=False))

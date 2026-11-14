from pathlib import Path
import importlib.util
import pandas as pd
import numpy as np

path=Path(__file__).resolve().parent/'run_all.py'
spec=importlib.util.spec_from_file_location('run_all_module', path)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
rows=[]
for code,cfg in m.CITY.items():
    loads={}
    for s,p in cfg['profiles'].items():
        d=m.read_load(p).rename(columns={'total_electric_kWh':'load_kWh'})
        b=m.monthly_bill_no_pv(d,cfg['distributor'],1)
        loads[s]={'hourly':d,'annual_kWh':d.load_kWh.sum(),'bill_RD':b.total_bill_RD.sum()}
    pv,_=m.pv_profile(cfg['epw']); gen=pv.pv_ac_kWh.to_numpy(); load=loads['E2']['hourly'].load_kWh.to_numpy()
    cap,_=m.knee_capacity(gen,load)
    # base E3/E4 bills only to complete loads dict; finance recalculates year-specific bills itself.
    for scen,ccap in [('E3',0),('E4',cap)]:
        r=m.dispatch(gen,load,ccap)
        h=loads['E2']['hourly'].copy(); h['pv_ac_kWh']=gen; h['grid_import_kWh']=r['import']; h['grid_export_kWh']=r['export']
        b=m.monthly_bill_net_meter(h,cfg['distributor'],cfg['ppa_usd_mwh'],1)
        loads[scen]={'bill_RD':b.total_bill_RD.sum()}
    base=loads['E0']['bill_RD']
    for scen in ['E1','E2','E3','E4']:
        f=m.finance(code,scen,base,loads,gen,cap)
        for year,cf in enumerate(f['cashflows_USD']):
            rows.append({
                'city':code,'name':cfg['name'],'scenario':scen,'year':year,
                'cashflow_USD':cf,
                'discounted_cashflow_USD':cf/((1+m.ASSUMPTIONS['discount_rate_real'])**year),
                'cumulative_cashflow_USD':sum(f['cashflows_USD'][:year+1]),
                'cumulative_discounted_cashflow_USD':sum(f['cashflows_USD'][i]/((1+m.ASSUMPTIONS['discount_rate_real'])**i) for i in range(year+1)),
                'bill_RD': (None if year==0 else f['year_bills_RD'][year-1]),
                'investment_USD':f['capex_USD'],'NPV_USD':f['NPV_USD'],'IRR_pct':f['IRR_pct']*1.0,
            })
out=Path('/mnt/data/Retrofit_RD_reproducible/outputs/summaries/financial_cashflows_25y.csv')
pd.DataFrame(rows).to_csv(out,index=False)
print(out, len(rows))

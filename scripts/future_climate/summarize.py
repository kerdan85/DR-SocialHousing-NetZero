"""Collect the present and future-climate results into fut/summary.json and fut/summary.csv.
Usage: python summarize.py <futroot>"""
import json, sys, csv
from pathlib import Path
FR = Path(sys.argv[1])
SC = ['present', 'ssp245_2050', 'ssp245_2080', 'ssp585_2050', 'ssp585_2080']
SITES = ['MC', 'SD', 'ST']; A = 2943.0; EF = 0.6277


def tmean(p):
    rows = [l.split(',') for l in open(p, encoding='latin-1').read().splitlines()[8:] if l.strip()]
    T = [float(r[6]) for r in rows]; Td = [float(r[7]) for r in rows]
    return sum(T) / len(T), sum(Td) / len(Td), sum(1 for t in T if t > 32)


out = {}; rows = []
for sc in SC:
    r = json.load(open(FR / 'runs' / sc / 'results/results.json'))
    out[sc] = {}
    for s in SITES:
        o = r[s]; F = o['full']; Nn = o['night']; fb = F['fin']['BTD']; nb = Nn['fin']['BTD']
        T, Td, h32 = tmean(FR / 'epw' / sc / f'{s}.epw')
        d = dict(scen=sc, site=s, T=T, Td=Td, h32=h32, PV=o['PV_MWh'],
                 cool0=F['E0']['cool_MWh'], cool1=F['E1']['cool_MWh'], cool2=F['E2']['cool_MWh'],
                 coolth0=F['E0']['cool_th_MWh'], lat0=F['E0']['latent_pct'], lat2=F['E2']['latent_pct'],
                 load0=F['E0']['load_MWh'], load2=F['E2']['load_MWh'],
                 nload0=Nn['E0']['load_MWh'], nload2=Nn['E2']['load_MWh'], ncool0=Nn['E0']['cool_MWh'], ncool2=Nn['E2']['cool_MWh'],
                 peak0=F['E0']['peak_cool_kW'], peak2=F['E2']['peak_cool_kW'],
                 eve0=F['E0']['evening_kW'], eve2=F['E2']['evening_kW'], eve4=F['E4']['evening_imp_kW'],
                 cap=F['E4']['cap'], SSR3=F['E3']['SSR'], SSR4=F['E4']['SSR'], imp3=F['E3']['imp_MWh'], imp4=F['E4']['imp_MWh'],
                 out48=F['E5']['ord_pct'], reserve=F['E5']['reserve'],
                 nSSR4=Nn['E4']['SSR'],
                 NPV_E2=fb['E2']['NPV'], NPV_E3=fb['E3']['NPV'], NPV_E4=fb['E4']['NPV'],
                 sav1_E2=fb['E2']['sav1'], nNPV_E2=nb['E2']['NPV'], nNPV_E3=nb['E3']['NPV'])
        d['red2'] = 100 * (1 - d['cool2'] / d['cool0'])
        d['ratio0'] = d['PV'] / d['load0']; d['ratio2'] = d['PV'] / d['load2']
        d['nratio0'] = d['PV'] / d['nload0']; d['nratio2'] = d['PV'] / d['nload2']
        d['impred'] = 100 * (1 - d['imp4'] / d['load0'])
        d['dPV_NPV'] = d['NPV_E3'] - d['NPV_E2']; d['dBat_NPV'] = d['NPV_E4'] - d['NPV_E3']
        d['ndPV_NPV'] = d['nNPV_E3'] - d['nNPV_E2']
        d['CO2_0'] = d['load0'] * EF; d['CO2_4'] = d['imp4'] * EF
        d['eui0'] = d['cool0'] * 1000 / A; d['eui2'] = d['cool2'] * 1000 / A
        out[sc][s] = d; rows.append(d)
json.dump(out, open(FR / 'summary.json', 'w'), indent=1)
with open(FR / 'summary.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
for d in rows:
    print(f"{d['scen']:12s} {d['site']} T{d['T']:5.2f} cool0 {d['cool0']:6.1f} cool2 {d['cool2']:6.1f} red {d['red2']:4.1f}% "
          f"PV {d['PV']:5.1f} r0 {d['ratio0']:.2f} r2 {d['ratio2']:.2f} nr0 {d['nratio0']:.2f} nr2 {d['nratio2']:.2f} "
          f"SSR4 {d['SSR4']:4.1f} cap {d['cap']:5.0f} imp-red {d['impred']:4.1f} out48 {d['out48']:5.1f} res {d['reserve']:4.0f} "
          f"NPV2 {d['NPV_E2']/1e3:6.0f} dPV {d['dPV_NPV']/1e3:5.0f} dBat {d['dBat_NPV']/1e3:5.0f} eve0 {d['eve0']:4.1f} eve2 {d['eve2']:4.1f}")

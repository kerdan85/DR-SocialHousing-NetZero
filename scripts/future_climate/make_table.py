"""Build the LaTeX table of future-climate results from fut/summary.json and the run results."""
import json, sys
from decimal import Decimal, ROUND_HALF_UP
FR=sys.argv[1]
S=json.load(open(f'{FR}/summary.json')); SC=['present','ssp245_2050','ssp585_2050','ssp245_2080','ssp585_2080']; X=['MC','SD','ST']
R={sc:json.load(open(f'{FR}/runs/{sc}/results/results.json')) for sc in SC}
def q(v,n=0): return str(Decimal(repr(v)).quantize(Decimal(1) if n==0 else Decimal('1.'+'0'*n),rounding=ROUND_HALF_UP))
def rg(vals,n=0):
    lo,hi=q(min(vals),n),q(max(vals),n)
    if lo==hi: return lo.replace('-','$-$')
    if min(vals)<0: return f'$-${lo[1:]} to $-${hi[1:]}' if max(vals)<0 else f'{lo} to {hi}'
    return f'{lo}--{hi}'
def row(lab,f,n=0,skip_present=False):
    cells=[]
    for sc in SC:
        if skip_present and sc=='present': cells.append('--'); continue
        cells.append(rg([f(sc,s) for s in X],n))
    return lab+' & '+' & '.join(cells)+r'\\'
P=S['present']
L=[]
L.append(row(r'Warming of mean air temperature (\degC)',lambda sc,s:S[sc][s]['T']-P[s]['T'],1,True))
L.append(row(r'Hours above 32\degC{} per year',lambda sc,s:S[sc][s]['h32']))
L.append(row(r'Cooling electricity, E0 (MWh/year)',lambda sc,s:S[sc][s]['cool0']))
L.append(row(r'Cooling electricity, E2 (MWh/year)',lambda sc,s:S[sc][s]['cool2']))
L.append(row(r'Cooling reduction by the retrofit (\%)',lambda sc,s:S[sc][s]['red2']))
L.append(row(r'PV/load ratio, E2, night-only',lambda sc,s:S[sc][s]['nratio2'],2))
L.append(row(r'PV/load ratio, E2, full conditioning',lambda sc,s:S[sc][s]['ratio2'],2))
L.append(row(r'Battery at the knee (kWh)',lambda sc,s:S[sc][s]['cap']))
L.append(row(r'Self-sufficiency ratio, E4 (\%)',lambda sc,s:S[sc][s]['SSR4']))
L.append(row(r'Grid imports avoided, E0 to E4 (\%)',lambda sc,s:S[sc][s]['impred']))
L.append(row(r'Evening load, E0 (kW)',lambda sc,s:S[sc][s]['eve0']))
L.append(row(r'Evening load, E2 (kW)',lambda sc,s:S[sc][s]['eve2']))
# reserve with critical load scaled
cells=[]
for sc in SC:
    v=[R[sc][s]['full']['E5']['reserve'] for s in X if R[sc][s]['full']['E5'].get('reserve_found',True)]
    nf=[s for s in X if not R[sc][s]['full']['E5'].get('reserve_found',True)]
    if not v: cells.append('n.s.')
    else: cells.append(rg(v)+('$^{a}$' if nf else ''))
L.append(r'SOC reserve for every 48~h outage, critical load 30\% of mean load (\%) & '+' & '.join(cells)+r'\\')
L.append(row(r'SOC reserve for every 48~h outage, present-day critical load (\%)',lambda sc,s:(R[sc][s]['full']['E5'].get('crit0_reserve') or R[sc][s]['full']['E5']['reserve'])))
L.append(row(r'NPV of the retrofit, E0 to E2 (thousand \usd)',lambda sc,s:S[sc][s]['NPV_E2']/1e3))
L.append(row(r'NPV of PV, E2 to E3 (thousand \usd)',lambda sc,s:S[sc][s]['dPV_NPV']/1e3))
L.append(row(r'NPV of storage, E3 to E4 (thousand \usd)',lambda sc,s:S[sc][s]['dBat_NPV']/1e3))
body='\n'.join(L)
tab=r'''\begin{table}[H]
\caption{Results under projected climates: ranges across the three sites, full conditioning and demand tariff unless stated. Warming is the change in mean annual air temperature relative to the TMYx files. The evening load is the mean of 18:00--23:00. NPV: net present value over 25 years; for PV and storage, the increment over the previous step. n.s.: no reserve serves every outage.}\label{tab:future}
\tablefont
\begin{widefloat}
\begin{tabularx}{\widewidth}{@{}L r r r r r@{}}
\toprule
 & & \multicolumn{2}{c}{\textbf{2050}} & \multicolumn{2}{c}{\textbf{2080}}\\
\cmidrule(lr){3-4}\cmidrule(l){5-6}
\textbf{Indicator} & \textbf{Present} & \textbf{SSP2-4.5} & \textbf{SSP5-8.5} & \textbf{SSP2-4.5} & \textbf{SSP5-8.5}\\
\midrule
'''+body+r'''
\bottomrule
\end{tabularx}
\end{widefloat}

\smallskip
{\footnotesize $^{a}$ Montecristi and Santiago; no reserve is sufficient in Santo Domingo.}
\end{table}
'''
open(f'{FR}/table_future.tex','w').write(tab); print(tab)

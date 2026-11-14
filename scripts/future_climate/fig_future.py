"""Figure: response of cooling, net-zero balance and outage reserve to warming (present + four CMIP6 futures).
Usage: python fig_future.py <futroot> <figdir> <nstyle_dir>"""
import sys, json
sys.path.insert(0, sys.argv[3])
import nstyle as N
import matplotlib.pyplot as plt
FR, OUT = sys.argv[1], sys.argv[2]
S = json.load(open(f'{FR}/summary.json'))
R = {sc: json.load(open(f'{FR}/runs/{sc}/results/results.json')) for sc in S}
SC = ['present', 'ssp245_2050', 'ssp585_2050', 'ssp245_2080', 'ssp585_2080']
MK = {'present': 'o', 'ssp245_2050': 's', 'ssp585_2050': 's', 'ssp245_2080': 'D', 'ssp585_2080': 'D'}
SITES = ['MC', 'SD', 'ST']

def dT(sc, s): return S[sc][s]['T'] - S['present'][s]['T']

f, ax = plt.subplots(1, 3, figsize=(N.W2, 2.25))
plt.subplots_adjust(wspace=0.48)
# a: cooling electricity
a = ax[0]
for s in SITES:
    x = [dT(sc, s) for sc in SC]
    a.plot(x, [S[sc][s]['cool0'] for sc in SC], '--', color=N.SITE[s], lw=0.9)
    a.plot(x, [S[sc][s]['cool0'] for sc in SC], 'o', ms=3.2, mfc='white', mec=N.SITE[s], mew=0.8)
    a.plot(x, [S[sc][s]['cool2'] for sc in SC], '-', color=N.SITE[s], lw=1.1)
    a.plot(x, [S[sc][s]['cool2'] for sc in SC], 'o', ms=3.2, color=N.SITE[s], mec='white', mew=0.4)
a.set_ylim(0, 480); a.set_xlim(-0.15, 3.75)
a.set_ylabel('Cooling electricity (MWh per year)')
a.set_xlabel('Warming of mean air temperature (°C)')
a.text(3.6, S['ssp585_2080']['MC']['cool0'] + 14, 'E0', fontsize=6.3, ha='right', color=N.INK2)
a.text(3.6, S['ssp585_2080']['MC']['cool2'] + 14, 'E2', fontsize=6.3, ha='right', color=N.INK2)
for s, y in zip(SITES, (0.95, 0.88, 0.81)):
    a.text(0.03, y, N.NAME[s], transform=a.transAxes, color=N.SITE[s], fontsize=6.3, fontweight='bold')
N.panel(a, 'a', x=-0.27)
a.set_title('Full conditioning', pad=6)

# b: annual PV / load
b = ax[1]
b.axhline(1.0, color='#1b7837', lw=0.8)
b.text(3.8, 1.0, 'net\nzero', fontsize=6.0, color='#1b7837', ha='left', va='center', clip_on=False, linespacing=1.0)
for s in SITES:
    x = [dT(sc, s) for sc in SC]
    b.plot(x, [S[sc][s]['nratio2'] for sc in SC], '-', color=N.SITE[s], lw=1.1)
    b.plot(x, [S[sc][s]['nratio2'] for sc in SC], 'o', ms=3.2, color=N.SITE[s], mec='white', mew=0.4)
    b.plot(x, [S[sc][s]['ratio2'] for sc in SC], ':', color=N.SITE[s], lw=1.0)
    b.plot(x, [S[sc][s]['ratio2'] for sc in SC], 'o', ms=2.6, mfc='white', mec=N.SITE[s], mew=0.7)
b.set_ylim(0.4, 1.55); b.set_xlim(-0.15, 3.75)
b.set_ylabel('Annual PV generation / load (E2)')
b.set_xlabel('Warming of mean air temperature (°C)')
b.text(0.1, 1.5, 'night-only', fontsize=6.0, color=N.INK2, va='center')
b.text(0.1, 0.47, 'full conditioning', fontsize=6.0, color=N.INK2, va='center')
N.panel(b, 'b', x=-0.27)
b.set_title('After the retrofit', pad=6)

# c: SOC reserve needed for every 48 h outage
c = ax[2]
for s in SITES:
    x = [dT(sc, s) for sc in SC]
    ys = [R[sc][s]['full']['E5']['reserve'] for sc in SC]
    ok = [R[sc][s]['full']['E5'].get('reserve_found', True) for sc in SC]
    y0 = [R[sc][s]['full']['E5'].get('crit0_reserve') or R[sc][s]['full']['E5']['reserve'] for sc in SC]
    c.plot(x, ys, '-', color=N.SITE[s], lw=1.1)
    for xi, yi, k in zip(x, ys, ok):
        if k: c.plot(xi, yi, 'o', ms=3.2, color=N.SITE[s], mec='white', mew=0.4)
        else: c.plot(xi, yi, 'x', ms=4.2, color=N.SITE[s], mew=1.0)
    c.plot(x, y0, '--', color=N.SITE[s], lw=0.9)
    c.plot(x, y0, 'o', ms=2.6, mfc='white', mec=N.SITE[s], mew=0.7)
c.set_ylim(25, 105); c.set_xlim(-0.15, 3.75)
c.set_ylabel('Reserve for every 48 h outage (% SOC)')
c.set_xlabel('Warming of mean air temperature (°C)')
c.text(3.7, 27, '—  critical load 30% of mean load\n- -  present-day critical load\n×   no reserve is sufficient', fontsize=5.6,
       ha='right', va='bottom', color=N.INK2, linespacing=1.35)
N.panel(c, 'c', x=-0.27)
c.set_title('Outage resilience (E5)', pad=6)
N.save(f, OUT, 'Figure_25_FutureClimate')
print('ok')

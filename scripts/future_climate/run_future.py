"""Future-climate runs. For each scenario folder fut/<scen>/ (built from the repository), run the six
EnergyPlus models per site that pipeline.py uses with the morphed weather file, then run pipeline.py.
Usage: python run_future.py <repo> <futroot> <scen> [<scen> ...]
Expects <futroot>/epw/<scen>/<SITE>.epw (SITE in MC, SD, ST)."""
import sys, os, shutil, subprocess
from pathlib import Path
REPO = Path(sys.argv[1]); FR = Path(sys.argv[2]); SCENS = sys.argv[3:]
EPWN = {'MC': 'DOM_MC_Monte.Cristi.784510_TMYx.epw', 'SD': 'DOM_SD_Las.Americas.Intl.AP.784850_TMYx.epw',
        'ST': 'DOM_ST_Santiago-Cibao.Intl.AP.784600_TMYx.epw'}
VP = REPO / 'inputs/energyplus/validation_package'; VAR = REPO / 'inputs/energyplus/variants_v25'
for sc in SCENS:
    D = FR / 'runs' / sc
    clim = D / 'inputs/energyplus/validation_package/climate'; prof = D / 'inputs/energyplus/hourly_profiles_v25'
    for p in (clim, prof, D / 'scripts', D / 'results'): p.mkdir(parents=True, exist_ok=True)
    for f in ('run_all.py', 'pipeline.py', 'analysis.py', 'run_energyplus.py'):
        shutil.copy(REPO / 'scripts' / f, D / 'scripts' / f)
    # Future climates: no state-of-charge reserve may give full 48 h coverage. Record the coverage with the
    # battery held full (reserve 100%) instead of failing; also report coverage for the present-day critical load.
    pp = D / 'scripts/pipeline.py'; t = pp.read_text()
    a = "crit=A['critical_load_fraction']*L2.mean(); n48=outage(gen,cap,r4['soc'],crit,48); best=minres(gen,L2,cap,crit)"
    b = ("crit=A['critical_load_fraction']*L2.mean(); n48=outage(gen,cap,r4['soc'],crit,48); best=minres(gen,L2,cap,crit)\n"
         "            found=best is not None\n"
         "            if not found: best=(1.0,)\n"
         "            full=outage(gen,cap,np.full(len(L2),cap),crit,48)\n"
         "            CRIT0=json.load(open(r'" + str(FR / 'runs/present/results/results.json') + "'))[code]['full']['E5']['crit_kW']\n"
         "            best0=minres(gen,L2,cap,CRIT0); full0=outage(gen,cap,np.full(len(L2),cap),CRIT0,48)")
    assert a in t; t = t.replace(a, b)
    a2 = "'fixed6kW_ord_pct':nf['success_pct'],"
    b2 = ("'fixed6kW_ord_pct':nf['success_pct'],'reserve_found':found,'full_soc_pct':full['success_pct'],'full_soc_served':full['mean_served_pct'],"
          "'crit0_kW':CRIT0,'crit0_reserve':(best0[0]*100 if best0 else None),'crit0_full_soc_pct':full0['success_pct'],"
          "'crit0_ord_pct':outage(gen,cap,r4['soc'],CRIT0,48)['success_pct'],")
    assert a2 in t; t = t.replace(a2, b2); pp.write_text(t)
    for c, n in EPWN.items():
        shutil.copy(FR / 'epw' / sc / f'{c}.epw', clim / n)
    jobs = []
    for c in EPWN:
        jobs += [(f'{c}_E0', VP / f'models/{c}_E0_final046.idf', 2.8), (f'{c}_E0_night', VAR / f'{c}_E0_night.idf', 2.8)]
        for e, cop in (('E1', 2.8), ('E2', 4.0)):
            for v in ('film', 'night'):
                jobs.append((f'{c}_{e}_{v}', VAR / f'{c}_{e}_{v}.idf', cop))
    # ground temperatures follow the change in annual mean air temperature of the morphed file
    def tmean(p):
        rows = [l.split(',') for l in open(p, encoding='latin-1').read().splitlines()[8:] if l.strip()]
        return sum(float(r[6]) for r in rows) / len(rows)
    dT = {c: tmean(clim / n) - tmean(FR / 'epw' / 'present' / f'{c}.epw') for c, n in EPWN.items()}
    (D / 'ground_dT.txt').write_text(repr(dT))
    def shift_ground(src, dst, d):
        lines = open(src, encoding='latin-1').read().split('\n'); out = []; inblk = False
        for l in lines:
            s = l.strip()
            if s.startswith('Site:GroundTemperature:BuildingSurface') or s.startswith('Site:GroundTemperature:FCfactorMethod'):
                inblk = True; out.append(l); continue
            if inblk and 'Ground Temperature' in l:
                val, rest = l.split('!-', 1); v = val.strip(); end = v[-1]; num = float(v[:-1])
                l = f'    {num + d:.2f}{end}' + ' ' * 20 + '!-' + rest
                if end == ';': inblk = False
            out.append(l)
        open(dst, 'w', encoding='latin-1').write('\n'.join(out))
    for name, idf, cop in jobs:
        out = prof / f'{name}_hourly.csv'
        if out.exists(): continue
        od = FR / 'work' / sc / name; od.mkdir(parents=True, exist_ok=True)
        gidf = od.parent / f'{name}_ground.idf'
        if sc == 'present': shutil.copy(idf, gidf)
        else: shift_ground(idf, gidf, dT[name[:2]])
        subprocess.run([sys.executable, str(REPO / 'scripts/run_energyplus.py'), str(gidf), str(clim / EPWN[name[:2]]), str(od), str(cop)], check=True)
        shutil.copy(od / 'hourly.csv', out); shutil.rmtree(od)
    subprocess.run([sys.executable, str(D / 'scripts/pipeline.py')], check=True, cwd=D)
    print(sc, 'pipeline done', flush=True)

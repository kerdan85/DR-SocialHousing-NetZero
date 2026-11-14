"""Run the 24 EnergyPlus simulations used in the paper (v25) and copy each hourly profile to
inputs/energyplus/hourly_profiles_v25/<run>_hourly.csv. COP = 2.8 for E0/E1, 4.0 for E2.
Usage: python scripts/run_simulations.py   (requires EnergyPlus 23.1.0; ~1-2 min per run)"""
from pathlib import Path
import subprocess, shutil, sys
ROOT = Path(__file__).resolve().parents[1]
VP = ROOT/'inputs/energyplus/validation_package'; VAR = ROOT/'inputs/energyplus/variants_v25'
OUT = ROOT/'inputs/energyplus/hourly_profiles_v25'; WORK = ROOT/'work'
EPW = {'MC':'DOM_MC_Monte.Cristi.784510_TMYx.epw','SD':'DOM_SD_Las.Americas.Intl.AP.784850_TMYx.epw','ST':'DOM_ST_Santiago-Cibao.Intl.AP.784600_TMYx.epw'}
jobs = []
for c in EPW:
    jobs.append((f'{c}_E0', VP/f'models/{c}_E0_final046.idf', 2.8))
    jobs.append((f'{c}_E0_night', VAR/f'{c}_E0_night.idf', 2.8))
    for e, cop in [('E1', 2.8), ('E2', 4.0)]:
        for v in ['film', 'VENT_film', 'night']:
            jobs.append((f'{c}_{e}_{v}', VAR/f'{c}_{e}_{v}.idf', cop))
for name, idf, cop in jobs:
    c = name[:2]; od = WORK/name
    subprocess.run([sys.executable, str(ROOT/'scripts/run_energyplus.py'), str(idf), str(VP/'climate'/EPW[c]), str(od), str(cop)], check=True)
    shutil.copy(od/'hourly.csv', OUT/f'{name}_hourly.csv')
print(len(jobs), 'runs done')

"""Rebuild the v25 IDF variants from the original models (validation_package/models).
  *_film      : E1/E2 glazing modelled as solar-control film on existing single glass (U = 5.8 W/m2K)
  *_VENT_film : same change applied to the natural-ventilation models
  *_night     : cooling availability restricted to 21:00-07:00 (E1/E2 night variants include the film change)
Usage: python scripts/make_variants.py"""
from pathlib import Path
import re
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT/'inputs/energyplus/validation_package/models'
DST = ROOT/'inputs/energyplus/variants_v25'; DST.mkdir(parents=True, exist_ok=True)
U_OLD = '    3,                        !- UFactor'
U_NEW = '    5.8,                      !- UFactor'
SCH = re.compile(r'(Cooling Availability Sch,.*?For: AllDays,[^\n]*\n)    Until: 24:00,             !- Field 3\n    1;                        !- Field 4', re.S)
NIGHT = r'\1    Until: 07:00, 1,\n    Until: 21:00, 0,\n    Until: 24:00, 1;'
def rd(p): return p.read_text(encoding='latin-1')
def wr(p, t): p.write_text(t, encoding='latin-1')
def film(t):
    assert t.count(U_OLD) == 1; return t.replace(U_OLD, U_NEW)
def night(t):
    t2, n = SCH.subn(NIGHT, t); assert n > 0; return t2
for c in ['MC', 'SD', 'ST']:
    for e in ['E1', 'E2']:
        wr(DST/f'{c}_{e}_film.idf', film(rd(SRC/f'{c}_{e}_final046.idf')))
        wr(DST/f'{c}_{e}_VENT_film.idf', film(rd(SRC/f'{c}_{e}_VENT_final046.idf')))
        wr(DST/f'{c}_{e}_night.idf', night(film(rd(SRC/f'{c}_{e}_final046.idf'))))
    wr(DST/f'{c}_E0_night.idf', night(rd(SRC/f'{c}_E0_final046.idf')))
print('variants written to', DST)

"""Build the future weather files used in the simulations from the Future Weather Generator (FWG) output.
Temperature, humidity, pressure, global horizontal irradiance, wind, sky cover and precipitation come from the
FWG morphed file. FWG recomputes the direct/diffuse split and the horizontal infrared radiation with its own
methods, which changes them even with a zero climate signal (checked with the FWG present-day control). To keep
only the climate signal: direct and diffuse (irradiance and illuminance) are the original values scaled by the
hourly ratio of morphed to original global irradiance, and infrared and opaque sky cover are the original values
plus the hourly difference between the morphed file and the FWG present-day control.
Usage: python build_epw.py <orig.epw> <fwg_future.epw> <fwg_control.epw> <out.epw>"""
import sys
o_p, f_p, c_p, out_p = sys.argv[1:5]
rd = lambda p: open(p, encoding='latin-1').read().splitlines()
O, F, C = rd(o_p), rd(f_p), rd(c_p)
hdr = F[:8]
ro = [l.split(',') for l in O[8:] if l.strip()]
rf = [l.split(',') for l in F[8:] if l.strip()]
rc = [l.split(',') for l in C[8:] if l.strip()]
assert len(ro) == len(rf) == len(rc) == 8760
out = []
for a, f, c in zip(ro, rf, rc):
    assert a[1:4] == f[1:4] == c[1:4], (a[:4], f[:4])
    r = list(f)
    g0, g1 = float(a[13]), float(f[13])
    k = g1 / g0 if g0 > 0 else 0.0
    for i in (14, 15, 16, 17, 18):          # DNI, DHI, global/direct/diffuse illuminance
        r[i] = str(int(round(float(a[i]) * k))) if i != 16 else str(int(round(float(a[16]) * k)))
    r[12] = str(int(round(float(a[12]) + float(f[12]) - float(c[12]))))                    # infrared
    r[23] = str(int(min(10, max(0, round(float(a[23]) + float(f[23]) - float(c[23]))))))   # opaque sky cover
    out.append(','.join(r))
note = 'COMMENTS 2,Direct/diffuse split and infrared from the original TMYx with the FWG climate signal applied (see build_epw.py)'
hdr = [h if not h.startswith('COMMENTS 2') else note for h in hdr]
open(out_p, 'w', encoding='latin-1', newline='\r\n').write('\n'.join(hdr + out) + '\n')

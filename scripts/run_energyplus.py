"""Run one EnergyPlus 23.1.0 simulation and write an hourly load profile in the same format as validation_package/hourly_profiles.
Usage: python scripts/run_energyplus.py <idf> <epw> <outdir> <COP>
Set ENERGYPLUS_DIR to the EnergyPlus install folder if energyplus is not on PATH."""
import sys,subprocess,os,re
import pandas as pd
idf,epw,out,cop=sys.argv[1],sys.argv[2],sys.argv[3],float(sys.argv[4])
EP=os.environ.get('ENERGYPLUS_DIR','')  # folder containing the energyplus executable (v23.1.0)
exe=os.path.join(EP,'energyplus') if EP else 'energyplus'
os.makedirs(out,exist_ok=True)
t=open(idf,encoding='latin-1').read()
add="""
Output:Variable,*,Zone Ideal Loads Supply Air Total Cooling Energy,hourly;
Output:Variable,*,Zone Ideal Loads Supply Air Sensible Cooling Energy,hourly;
Output:Variable,*,Zone Ideal Loads Supply Air Latent Cooling Energy,hourly;
Output:Meter,InteriorLights:Electricity,hourly;
Output:Meter,InteriorEquipment:Electricity,hourly;
"""
open(out+'/in.idf','w',encoding='latin-1').write(t+add)
r=subprocess.run([exe,'-r','-w',epw,'-d',out,out+'/in.idf'],capture_output=True,text=True)
if r.returncode: print(r.stdout[-2000:],r.stderr[-2000:]); sys.exit(1)
d=pd.read_csv(out+'/eplusout.csv')
J=3.6e6
tot=d[[c for c in d.columns if 'Supply Air Total Cooling Energy' in c]].sum(axis=1)/J
sen=d[[c for c in d.columns if 'Supply Air Sensible Cooling Energy' in c]].sum(axis=1)/J
lat=d[[c for c in d.columns if 'Supply Air Latent Cooling Energy' in c]].sum(axis=1)/J
li=d[[c for c in d.columns if c.startswith('InteriorLights:Electricity') and 'Hourly' in c]].iloc[:,0]/J
eq=d[[c for c in d.columns if c.startswith('InteriorEquipment:Electricity') and 'Hourly' in c]].iloc[:,0]/J
dt=d['Date/Time'].str.strip()
mo=dt.str[0:2].astype(int); dy=dt.str[3:5].astype(int); hr=dt.str[7:9].astype(int)
p=pd.DataFrame({'hour_index':range(1,len(d)+1),'year':2002,'month':mo,'day':dy,'hour_ending':hr,'minute':0,'day_type':'',
 'cooling_thermal_kWh':tot,'cooling_electric_kWh':tot/cop,'lighting_kWh':li,'equipment_kWh':eq})
p['total_electric_kWh']=p.cooling_electric_kWh+p.lighting_kWh+p.equipment_kWh
p['sensible_kWh']=sen; p['latent_kWh']=lat
p.to_csv(out+'/hourly.csv',index=False)
print(os.path.basename(idf), 'thermal',round(tot.sum()),'sens',round(sen.sum()),'lat',round(lat.sum()),'elec',round(p.cooling_electric_kWh.sum()),'peak',round(p.cooling_electric_kWh.max(),1),'lights',round(li.sum()),'eq',round(eq.sum()))

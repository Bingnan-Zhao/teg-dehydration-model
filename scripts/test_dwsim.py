"""
test_dwsim.py  -  check that Python can drive DWSIM.

Mixes 100 kg/s of water at 300 K with 50 kg/s at 348 K.
Expected outlet ~316 K. Result is printed and saved to results/test_dwsim.txt.
Run:  python test_dwsim.py
"""
import os
from dwsim_setup import load_dwsim, PROJECT

Flowsheet, PropertyPackages, Q = load_dwsim()

fs = (Flowsheet.Create("test_mixer")
      .WithCompound("Water")
      .WithPropertyPackage(PropertyPackages.SteamTables))
a = fs.AddMaterialStream("cold").At(Q.Kelvin(300.0), Q.Pascal(101325.0)).WithMassFlow(Q.KgPerSecond(100.0))
b = fs.AddMaterialStream("hot").At(Q.Kelvin(348.0), Q.Pascal(101325.0)).WithMassFlow(Q.KgPerSecond(50.0))
out = fs.AddMaterialStream("out")
fs.AddMixer("MIX-1").ConnectFeed(a, 0).ConnectFeed(b, 1).ConnectProduct(out, 0)
fs.Solve()

T = out.TemperatureK
ok = abs(T - 316.0) < 1.0
msg = f"Outlet T = {T:.2f} K ({T - 273.15:.2f} C), flow = {out.MassFlowKgPerSecond:.1f} kg/s -> {'PASS' if ok else 'CHECK'}"
print(msg)
os.makedirs(os.path.join(PROJECT, "results"), exist_ok=True)
open(os.path.join(PROJECT, "results", "test_dwsim.txt"), "w").write(msg + "\n")

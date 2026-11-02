"""
dwsim_setup.py  -  connect Python to an installed DWSIM 10 (Fluent API).

Usage in other scripts:
    from dwsim_setup import load_dwsim
    Flowsheet, PropertyPackages, Q = load_dwsim()

DWSIM 10 ships its own .NET 10 runtime inside its install folder. pythonnet
needs that runtime laid out like a normal .NET install, so on first run this
script copies the runtime files into  <project>/.dnroot  (one time, ~100 MB).
Set the environment variable DWSIM_BIN to the DWSIM folder if it is not found
automatically.
"""
import glob
import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
DNROOT = os.path.join(PROJECT, ".dnroot")
WIN = sys.platform.startswith("win")


def find_dwsim_bin():
    env = os.environ.get("DWSIM_BIN")
    if env and os.path.exists(os.path.join(env, "DWSIM.Automation.FluentAPI.dll")):
        return env
    roots = []
    if WIN:
        for var in ("LOCALAPPDATA", "ProgramFiles", "ProgramFiles(x86)", "APPDATA"):
            v = os.environ.get(var)
            if v:
                roots += [v, os.path.join(v, "Programs")]
    else:
        roots += ["/usr/lib", "/opt", os.path.expanduser("~")]
    for r in roots:
        for pat in ("DWSIM*", "DWSIM*/*", "DWSIM*/*/*"):
            for d in glob.glob(os.path.join(r, pat)):
                if os.path.exists(os.path.join(d, "DWSIM.Automation.FluentAPI.dll")):
                    return d
    raise FileNotFoundError(
        "DWSIM folder not found. Set DWSIM_BIN to the folder that contains "
        "DWSIM.Automation.FluentAPI.dll")


def build_dnroot(bin_dir):
    """Lay out DWSIM's bundled runtime as <root>/shared/Microsoft.NETCore.App/<ver>."""
    deps = None
    for f in glob.glob(os.path.join(bin_dir, "*.deps.json")):
        d = json.load(open(f, encoding="utf-8-sig"))
        tname = d["runtimeTarget"]["name"]
        for key, entry in d["targets"][tname].items():
            if key.startswith("runtimepack.Microsoft.NETCore.App.Runtime"):
                deps = (tname, key, entry)
                break
        if deps:
            break
    if not deps:
        raise RuntimeError("No .NET runtime pack found in the DWSIM folder.")
    tname, key, entry = deps
    ver = key.split("/")[-1]
    shared = os.path.join(DNROOT, "shared", "Microsoft.NETCore.App", ver)
    marker = os.path.join(shared, "Microsoft.NETCore.App.deps.json")
    if not os.path.exists(marker):
        print(f"First run: copying .NET {ver} runtime from DWSIM (one time)...")
        os.makedirs(shared, exist_ok=True)
        for sect in ("runtime", "native"):
            for name in entry.get(sect, {}):
                src = os.path.join(bin_dir, name)
                if os.path.exists(src):
                    shutil.copy2(src, os.path.join(shared, os.path.basename(name)))
        fxr = "hostfxr.dll" if WIN else "libhostfxr.so"
        fxr_dir = os.path.join(DNROOT, "host", "fxr", ver)
        os.makedirs(fxr_dir, exist_ok=True)
        shutil.copy2(os.path.join(bin_dir, fxr), fxr_dir)
        pack = key.replace("runtimepack.", "")
        json.dump({"runtimeTarget": {"name": tname, "signature": ""},
                   "compilationOptions": {},
                   "targets": {tname: {pack: entry}},
                   "libraries": {pack: {"type": "package", "serviceable": True, "sha512": ""}}},
                  open(marker, "w"))
        json.dump({"runtimeOptions": {"tfm": "net10.0", "configProperties": {}}},
                  open(os.path.join(shared, "Microsoft.NETCore.App.runtimeconfig.json"), "w"))
    rc = os.path.join(DNROOT, "app.runtimeconfig.json")
    json.dump({"runtimeOptions": {"tfm": "net10.0",
                                  "framework": {"name": "Microsoft.NETCore.App", "version": ver},
                                  "configProperties": {
                                      "System.Runtime.Serialization.EnableUnsafeBinaryFormatterSerialization": True}}},
              open(rc, "w"))
    return rc


def load_dwsim():
    bin_dir = find_dwsim_bin()
    rc = build_dnroot(bin_dir)
    from pythonnet import load
    load("coreclr", runtime_config=rc, dotnet_root=DNROOT)
    import clr
    sys.path.append(bin_dir)
    clr.AddReference("DWSIM.Automation.FluentAPI")
    from DWSIM.Automation.FluentAPI import Flowsheet, PropertyPackages, Q
    return Flowsheet, PropertyPackages, Q

#!/usr/bin/env python3
"""Regression suite: short NHFOAM runs whose gauge records are compared with
stored references (regression_ref.json), plus parallel and restart checks.

usage: regression.py [--update] [--np N] [--keep] <workDir> [case ...]

    --update   store the results as the new references
    --np N     processors of the parallel check (default 2; 0: skip it)

Each case writes a fingerprint: for up to 8 gauges the elevation at the end
and its maximum and minimum over the run, the shoreline (wetting and drying
runs) and the water volume change of closed runs.  A case passes when every
number is within its absolute tolerance of the reference.  On one machine
the runs repeat to round-off; the tolerances allow for other compilers and
platforms, and are wide for the breaking cases, where the on/off breaking
flags turn a perturbation of 1e-10 into millimetres within seconds.  Two checks need
no reference: the parallel run (decomposePar, mpirun) and the restart (run
to half time, restart from the written fields) must reproduce the serial,
continuous run to round-off.  The OpenFOAM environment must be sourced and
NHFOAM built; gmsh is needed for the unstructured case (skipped without it).
Run time: about 3 minutes on one core.
"""
import json, os, re, shutil, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_case as mc
import wetdry

REF = os.path.join(HERE, "regression_ref.json")


# ------------------------------------------------------------------ cases
def c_standing(case):
    mc.standing_case(case, kh=1.0, nLayers=3, Nx=32, nPeriods=2, stepsPerPeriod=100)


def c_solitary(case):
    mc.solitary_case(case, H=0.1, h=1.0, nLayers=5, dx=0.1, dt=0.02, L=30.0, x0=5.0,
                     endTime=4.0)


def c_bar_bc(case):
    # gauges near the wavemaker: the waves reach about 8 m in 10 s
    mc.bar_case(case, nLayers=5, dx=0.04, dt=0.01, endTime=10.0, H=0.041, T=1.01,
                wavemaker="bc", returnCurrent=True, gauges=[1.5, 3.0, 4.5])


def c_irregular(case):
    import wavemaker as wm
    wm.DX, wm.NZ = 0.04, 5
    wm.flume(case, 10.0, 8.0, [1.0, 3.0, 5.0], zones=[wm.zone("outlet", 7, 10)],
             left=wm.bc("left", wm.IRR, ramp=1.5))


def c_dambreak(case):
    wetdry.dambreak(case, dx=0.05, nz=3, tend=1.0)


def c_runup(case):
    wetdry.runup_nhfoam(case, Cd=0.0, tstar_end=30.0)


def c_tk(case):
    wetdry.tingkirby(case, nz=5, dx=0.05, dt=0.01, nPeriods=6, brk_model="hybrid",
                     alpha=0.6, beta=0.3,
                     brk_extra="\n    betaFinal   0.15;\n    Tstar       5;\n    delta       1.2;")


def c_kw(case):
    wetdry.tingkirby(case, nz=10, dx=0.05, dt=0.01, nPeriods=1, breaking=False,
                     rans=dict(model="kOmegaLF"))


def c_uplift(case):
    mc.uplift_case(case)
    cd = os.path.join(case, "system/controlDict")
    txt = open(cd).read()
    end = float(re.search(r"\nendTime\s+([^;]+);", txt).group(1))
    open(cd, "w").write(re.sub(r"\nendTime\s+[^;]+;", f"\nendTime         {0.25 * end:.6g};", txt))


def c_strip(case):
    import unstructured
    unstructured.runup_strip(case, W=0.08, lc=0.04, tstar_end=20.0, dt=0.004)


def c_ts_eta(case):
    # solitary wave (H/h = 0.3) from the record of its elevation at the inflow
    import timeseries_input as ti
    ti.flat_case(case, "eta", Hh=0.3, dx=0.04, nz=3, dt=0.01, endTime=5.0)


def c_ts_paddle(case):
    # the same wave from a piston-paddle record, without absorption
    import timeseries_input as ti
    ti.flat_case(case, "paddle", Hh=0.3, dx=0.04, nz=3, dt=0.01, endTime=5.0)


CASES = {
    "standing": (c_standing, 1e-8),
    "solitary": (c_solitary, 1e-7),
    "bar_bc": (c_bar_bc, 1e-6),
    "irregular": (c_irregular, 1e-6),
    "dambreak": (c_dambreak, 1e-6),
    "runup": (c_runup, 1e-3),
    "tk": (c_tk, 5e-3),
    "kw": (c_kw, 1e-5),
    "uplift": (c_uplift, 1e-8),
    "strip": (c_strip, 1e-4),
    "ts_eta": (c_ts_eta, 1e-6),
    "ts_paddle": (c_ts_paddle, 1e-6),
}


# ------------------------------------------------------------------ running
def run(case, cmd, log):
    with open(os.path.join(case, log), "w") as f:
        r = subprocess.run(cmd, cwd=case, stdout=f, stderr=subprocess.STDOUT)
    return r.returncode == 0


def run_case(case):
    """Mesh, initial conditions and NHFOAM, as the case's Allrun would."""
    allrun = open(os.path.join(case, "Allrun")).read() if os.path.exists(
        os.path.join(case, "Allrun")) else ""
    if "blockMesh" in allrun and not os.path.exists(os.path.join(case, "constant/polyMesh/points")):
        if not run(case, ["blockMesh"], "log.blockMesh"):
            return False
    if "nhSetWaves" in allrun and not run(case, ["nhSetWaves"], "log.nhSetWaves"):
        return False
    return run(case, ["NHFOAM"], "log.NHFOAM")


def set_dict(case, rel, key, value):
    p = os.path.join(case, rel)
    txt = open(p).read()
    new, n = re.subn(rf"\n{key}\s+[^;]+;", f"\n{key}    {value};", txt)
    if n == 0:
        new = txt.replace("\napplication", f"\n{key}    {value};\n\napplication", 1)
    open(p, "w").write(new)


# ------------------------------------------------------------------ fingerprints
def gauge_file(case, start="0"):
    return os.path.join(case, "postProcessing/gauges", start, "eta.dat")


def fingerprint(case):
    out = {}
    gf = gauge_file(case)
    if os.path.exists(gf):
        d = np.loadtxt(gf, comments="#", ndmin=2)
        n = d.shape[1] - 1
        idx = np.unique(np.linspace(0, n - 1, min(n, 8)).round().astype(int))
        for i in idx:
            e = d[:, i + 1]
            out[f"g{i}_end"] = float(e[-1])
            out[f"g{i}_max"] = float(e.max())
            out[f"g{i}_min"] = float(e.min())
    sf = os.path.join(case, "postProcessing/wettingDrying/0/shoreline.dat")
    if os.path.exists(sf):
        s = np.loadtxt(sf, comments="#", ndmin=2)
        out["shore_max"] = float(s[:, 2].max())
        out["shore_end"] = float(s[-1, 2])
    log = open(os.path.join(case, "log.NHFOAM")).read()
    v = re.findall(r"water volume change = ([-0-9.eE+]+)", log)
    if v:
        out["volume_end"] = float(v[-1])
    return out


def internal(case, time, name):
    """Internal field values of a scalar field."""
    txt = open(os.path.join(case, time, name)).read()
    i = txt.index("internalField")
    m = re.search(r"nonuniform\s+List<scalar>\s*(\d+)\s*\(", txt[i:])
    a = i + m.end()
    return np.array(txt[a:txt.index(")", a)].split(), float)


def last_gauges(case, start):
    d = np.loadtxt(gauge_file(case, start), comments="#", ndmin=2)
    return d[-1, 0], d[-1, 1:]


# ------------------------------------------------------------------ checks
def check_parallel(work, base, nproc):
    """The tk case (breaking, wetting and drying) on nproc processors."""
    case = os.path.join(work, "tk_par")
    shutil.rmtree(case, ignore_errors=True)
    c_tk(case)
    mc.write(case, "system/decomposeParDict", "dictionary", f"""numberOfSubdomains {nproc};
method          structured;     // whole columns on one processor
structuredCoeffs
{{
    method      simple;
    patches     (bottom);
    simpleCoeffs
    {{
        n       ({nproc} 1 1);
    }}
}}
""")
    if not run(case, ["blockMesh"], "log.blockMesh"):
        return False, "blockMesh failed"
    if not run(case, ["decomposePar", "-force"], "log.decomposePar"):
        return False, "decomposePar failed"
    t0 = time.time()
    ok = run(case, ["mpirun", "--allow-run-as-root", "--oversubscribe", "-np", str(nproc),
                    "NHFOAM", "-parallel"], "log.NHFOAM")
    wall = time.time() - t0
    if not ok:
        return False, "parallel run failed"
    ts, es = last_gauges(base, "0")
    tp, ep = last_gauges(case, "0")
    diff = float(np.abs(es - ep).max())
    ok = abs(ts - tp) < 1e-9 and diff < 1e-6
    return ok, f"max |eta_serial - eta_parallel| at t = {tp:g} s: {diff:.2e} m ({wall:.0f} s wall)"


def check_restart(work, base):
    """The tk case run to half time, then restarted from the written fields."""
    case = os.path.join(work, "tk_rst")
    shutil.rmtree(case, ignore_errors=True)
    c_tk(case)
    cd = os.path.join(case, "system/controlDict")
    txt = open(cd).read()
    end = float(re.search(r"\nendTime\s+([^;]+);", txt).group(1))
    dt = float(re.search(r"\ndeltaT\s+([^;]+);", txt).group(1))
    half = round(0.5 * end / dt) * dt
    # the mesh as in the continuous run (blockMesh writes the points with
    # writePrecision: a different precision is a different mesh)
    if not run(case, ["blockMesh"], "log.blockMesh"):
        return False, "blockMesh failed"
    set_dict(case, "system/controlDict", "endTime", f"{half:.10g}")
    set_dict(case, "system/controlDict", "writeInterval", str(int(round(half / dt))))
    # exact restarts need the fields at full precision (17 digits, or binary)
    set_dict(case, "system/controlDict", "writePrecision", "17")
    if not run_case(case):
        return False, "first half failed"
    set_dict(case, "system/controlDict", "startFrom", "latestTime")
    set_dict(case, "system/controlDict", "endTime", f"{end:.10g}")
    if not run(case, ["NHFOAM"], "log.NHFOAM.restart"):
        return False, "restart failed"
    starts = sorted(os.listdir(os.path.join(case, "postProcessing/gauges")), key=float)
    ts, es = last_gauges(base, "0")
    tr, er = last_gauges(case, starts[-1])
    diff = float(np.abs(es - er).max())
    tname = f"{end:g}"
    dfield = float(np.abs(internal(base, tname, "eta") - internal(case, tname, "eta")).max())
    ok = abs(ts - tr) < 1e-9 and diff < 1e-9 and dfield < 1e-9
    return ok, (f"max |eta_continuous - eta_restart| at t = {tr:g} s: gauges {diff:.2e} m,"
                f" field {dfield:.2e} m")


# ------------------------------------------------------------------ main
def main(argv):
    update = "--update" in argv
    keep = "--keep" in argv
    nproc = 2
    if "--np" in argv:
        i = argv.index("--np")
        nproc = int(argv[i + 1])
        del argv[i:i + 2]
    argv = [a for a in argv if a not in ("--update", "--keep")]
    work = os.path.abspath(argv[0])
    names = argv[1:] or list(CASES)
    os.makedirs(work, exist_ok=True)
    ref = json.load(open(REF)) if os.path.exists(REF) else {}
    have_gmsh = shutil.which("gmsh") is not None
    results, failed = {}, []
    for n in names:
        make, atol = CASES[n]
        if n == "strip" and not have_gmsh:
            print(f"{n:10s} SKIP (no gmsh)")
            continue
        case = os.path.join(work, n)
        shutil.rmtree(case, ignore_errors=True)
        t0 = time.time()
        make(case)
        ok = run_case(case)
        wall = time.time() - t0
        if not ok:
            print(f"{n:10s} FAIL (run failed, see {case})")
            failed.append(n)
            continue
        fp = fingerprint(case)
        results[n] = dict(atol=atol, values=fp)
        if n in ref and not update:
            rv = ref[n]["values"]
            worst, wk = 0.0, ""
            for k, v in rv.items():
                dv = abs(fp.get(k, np.inf) - v)
                if dv > worst:
                    worst, wk = dv, k
            ok = worst <= ref[n]["atol"]
            what = f"largest difference {worst:.2e} ({wk})" if worst > 0 else "identical"
            print(f"{n:10s} {'PASS' if ok else 'FAIL'}  {what}, tolerance"
                  f" {ref[n]['atol']:.0e}  {wall:.0f} s")
            if not ok:
                failed.append(n)
        else:
            print(f"{n:10s} {'stored' if update else 'no reference'}  {wall:.0f} s")
    if "tk" in names and os.path.exists(os.path.join(work, "tk", "log.NHFOAM")):
        for label, (ok, msg) in (
            ("restart", check_restart(work, os.path.join(work, "tk"))),
            *((("parallel", check_parallel(work, os.path.join(work, "tk"), nproc)),)
              if nproc > 1 else ())):
            print(f"{label:10s} {'PASS' if ok else 'FAIL'}  {msg}")
            if not ok:
                failed.append(label)
    if update:
        ref.update(results)
        with open(REF, "w") as f:
            json.dump(ref, f, indent=1, sort_keys=True)
        print("references written to", REF)
    if not keep and not failed:
        shutil.rmtree(work, ignore_errors=True)
    print("all passed" if not failed else "failed: " + " ".join(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

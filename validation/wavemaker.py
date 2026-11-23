#!/usr/bin/env python3
"""Wave-input tests in a flat flume (h = 0.4 m, 10 layers, dx = 2 cm):

  reg_bc, reg_zone      regular waves (H = 2 cm, T = 1.5 s) from the boundary
                        wavemaker nhWaveVelocity and from a relaxation zone,
                        sponge at the far end
  wall_abs, wall_noabs, wall_noabs_norc
                        boundary wavemaker in a closed 10 m flume (wall at the
                        far end): with and without active absorption, with and
                        without the return current
  zone_wall_rc, zone_wall_norc
                        relaxation zone (0-6 m) in a closed flume, wall at
                        16 m: target with and without the return current
  abs_T<T>_<celerity>, sponge_T<T>
                        absorbing boundary (nhWaveVelocity with zero wave
                        height) at the far end, celerity sqrt(gh) (shallow) or
                        the phase speed (wave), against a sponge layer, for
                        T = 1.0, 1.5, 2.5 s (kh = 1.72, 0.96, 0.53)
  irr_zone, irr_bc      JONSWAP waves (Hs = 3 cm, Tp = 1.5 s, gamma = 3.3) from
                        a relaxation zone and from the boundary wavemaker
                        (TVD convection, limitedLinearV)
  irr_bc_flt            irr_bc with filtered central convection

usage: wavemaker.py <root> [case ...]     (all cases by default)
"""
import os, sys
from make_case import (constant_files, system_files, blockmesh_flat, common_fields,
                       wave_bc, set_patch)

h, H = 0.4, 0.02
DX, DT, NZ = 0.02, 0.01, 10
IRR = """theory      irregular;
        Hs          0.03;
        Tp          1.5;
        gamma       3.3;
        nFrequencies 100;
        seed        1;"""


def regular(T, H=H):
    return f"""theory      linear;
        waveHeight  {H:g};
        wavePeriod  {T:g};"""


def wave_dict(theory, ramp, extra="", indent=8):
    pad = " " * indent
    body = f"""{theory}
depth       {h:g};
direction   (1 0 0);
rampTime    {ramp:g};
{extra}"""
    lines = [l.strip() for l in body.strip().splitlines() if l.strip()]
    return "\n".join(pad + l for l in lines)


def bc(patch, theory, ramp=3.0, absorb=True, celerity=None, extra=""):
    return wave_bc(patch, wave_dict(theory, ramp, extra), absorb=absorb, celerity=celerity)


def zone(name, start, end, theory=None, ramp=3.0, extra=""):
    if theory is None:
        body = "        type        sponge;\n"
    else:
        body = "        type        wave;\n" + wave_dict(theory, ramp, extra) + "\n"
    return f"""    {name}
    {{
{body}        start       ({start:g} 0 0);
        end         ({end:g} 0 0);
        sigmaMax    20;
        exponent    2;
    }}
"""


def flume(case, L, endTime, gauges, zones=(), left=None, right=None, conv=None):
    Nx = int(round(L / DX))
    gtxt = "\n".join(f"    G{i:03d}    ({x:g} {DX/2:g} 0);" for i, x in enumerate(gauges))
    ztxt = ""
    if zones:
        ztxt = "\nforcingZones\n{\n" + "".join(zones) + "}\n"
    constant_files(case, nh_extra=f"\ngauges\n{{\n{gtxt}\n}}\n{ztxt}")
    nsteps = int(round(endTime / DT))
    system_files(case, endTime=endTime, deltaT=DT, writeInterval=nsteps, nCorr=3,
                 pTol=1e-10, libs=True, conv=conv)
    blockmesh_flat(case, L, h, DX, Nx, NZ)
    common_fields(case)
    if left:
        set_patch(case, "U", "left", left)
    if right:
        set_patch(case, "U", "right", right)
    with open(os.path.join(case, "Allrun"), "w") as f:
        f.write("#!/bin/sh\ncd \"${0%/*}\" || exit 1\nblockMesh > log.blockMesh 2>&1\n"
                "NHFOAM > log.NHFOAM 2>&1\n")
    os.chmod(os.path.join(case, "Allrun"), 0o755)


def frange(a, b, d):
    n = int(round((b - a) / d))
    return [round(a + i * d, 6) for i in range(n + 1)]


def cases(root):
    c = {}
    g = frange(0.5, 13.5, 0.25)
    c["reg_bc"] = lambda p: flume(p, 20.0, 60.0, g, zones=[zone("outlet", 14, 20)],
                                  left=bc("left", regular(1.5)))
    c["reg_zone"] = lambda p: flume(p, 26.0, 60.0, [6 + x for x in g],
                                    zones=[zone("inlet", 6, 0, regular(1.5)),
                                           zone("outlet", 20, 26)])
    gw = frange(0.25, 9.75, 0.25)
    c["wall_abs"] = lambda p: flume(p, 10.0, 90.0, gw, left=bc("left", regular(1.5)))
    c["wall_noabs"] = lambda p: flume(p, 10.0, 90.0, gw,
                                      left=bc("left", regular(1.5), absorb=False))
    c["wall_noabs_norc"] = lambda p: flume(p, 10.0, 90.0, gw,
                                           left=bc("left", regular(1.5), absorb=False,
                                                   extra="returnCurrent no;"))
    gz = frange(6.25, 15.75, 0.25)
    for tag, rc in (("rc", ""), ("norc", "returnCurrent no;")):
        c[f"zone_wall_{tag}"] = (lambda p, rc=rc:
            flume(p, 16.0, 90.0, gz, zones=[zone("inlet", 6, 0, regular(1.5), extra=rc)]))
    ga = frange(6.5, 13.5, 0.1)
    for T in (1.0, 1.5, 2.5):
        tag = f"{T:g}".replace(".", "p")
        gen = zone("inlet", 6, 0, regular(T), ramp=2 * T)
        for cel in ("shallow", "wave"):
            c[f"abs_T{tag}_{cel}"] = (lambda p, T=T, gen=gen, cel=cel:
                flume(p, 14.0, 40 * T, ga, zones=[gen],
                      right=bc("right", regular(T, H=0.0), celerity=cel)))
        c[f"sponge_T{tag}"] = (lambda p, T=T, gen=gen:
            flume(p, 22.0, 40 * T, ga, zones=[gen, zone("outlet", 14, 22)]))
    gi = [8.0, 10.0, 12.0, 14.0, 16.0, 18.0]
    c["irr_zone"] = lambda p: flume(p, 26.0, 230.0, gi,
                                    zones=[zone("inlet", 6, 0, IRR),
                                           zone("outlet", 20, 26)])
    c["irr_bc"] = lambda p: flume(p, 20.0, 230.0, [x - 6 for x in gi],
                                  zones=[zone("outlet", 14, 20)],
                                  left=bc("left", IRR))
    c["irr_bc_flt"] = lambda p: flume(p, 20.0, 230.0, [x - 6 for x in gi],
                                      zones=[zone("outlet", 14, 20)],
                                      left=bc("left", IRR), conv="Gauss filteredLinear2V 0.2 0")
    return c


if __name__ == "__main__":
    root = sys.argv[1]
    allc = cases(root)
    names = sys.argv[2:] or list(allc)
    for n in names:
        allc[n](os.path.join(root, n))
        print("wrote", os.path.join(root, n))

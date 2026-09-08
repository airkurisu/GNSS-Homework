"""Compare broadcast states and SPP positions with user-supplied SP3/CLK/SINEX.

No third-party Python packages. Products stay outside the repository.
Orbit comparison uses exact SP3 epochs, never orbit interpolation/extrapolation.
Clock polynomials (without the periodic Kepler relativistic term) are compared
with raw CLK corrections; raw and epoch/system-demeaned differences are retained.
"""
import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import math
from pathlib import Path
import statistics
import subprocess
from collections import defaultdict

C = 299792458.0
OMEGA = 7.2921151467e-5
GPS_ORIGIN = dt.datetime(1980, 1, 6)


def gpstime(fields):
    y, m, d, h, minute = map(int, fields[:5])
    return round((dt.datetime(y, m, d, h, minute)-GPS_ORIGIN).total_seconds()+float(fields[5]), 6)


def norm(v):
    return math.sqrt(sum(x*x for x in v))


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def unit(v):
    length = norm(v)
    return [x/length for x in v]


def stats(values):
    values = [x for x in values if math.isfinite(x)]
    if not values:
        return {"count": 0}
    return {"count": len(values), "mean": statistics.mean(values),
            "rms": math.sqrt(statistics.mean(x*x for x in values)),
            "max_abs": max(map(abs, values))}


def load_sp3(path):
    records, times = {}, []
    header = []
    time = None
    with open(path, encoding="ascii") as source:
        for line in source:
            if line.startswith("*"):
                time = gpstime(line[1:].split())
                times.append(time)
            elif time is None:
                header.append(line.rstrip())
            elif line.startswith("P"):
                sat = line[1:4]
                if sat[0] not in "GCER":
                    continue
                xyz = [float(line[i:i+14])*1000 for i in (4, 18, 32)]
                if norm(xyz) == 0 or any(abs(x) >= 999999000 for x in xyz):
                    continue
                if len(line) > 79 and line[79] == "P":
                    continue  # predicted orbit
                records[(sat, time)] = xyz
    if not any(line.startswith("%c") and "GPS" in line for line in header):
        raise ValueError("SP3 must explicitly use GPS time")
    if len(times) < 9 or len(set(times)) != len(times):
        raise ValueError("Need at least nine distinct SP3 epochs")
    return records, sorted(times), header


def load_clk(path):
    result, header = {}, []
    in_header = True
    with open(path, encoding="ascii") as source:
        for line in source:
            if in_header:
                header.append(line.rstrip())
                if "END OF HEADER" in line:
                    in_header = False
                continue
            if not line.startswith("AS "):
                continue
            fields = line.split()
            sat, time = fields[1], gpstime(fields[2:8])
            clock = float(fields[9].replace("D", "E"))
            if sat[0] in "GCER" and math.isfinite(clock):
                result[(sat, time)] = clock
    if not any("TIME SYSTEM ID" in line and "GPS" in line[:10] for line in header):
        raise ValueError("CLK must explicitly use GPS time")
    return result, header


def load_sinex(path, station):
    axes, epochs, uncertainties = {}, set(), {}
    block = ""
    with open(path, encoding="ascii") as source:
        for line in source:
            if line.startswith("+"):
                block = line.strip()[1:]
            elif line.startswith("-"):
                block = ""
            elif block == "SOLUTION/ESTIMATE" and not line.startswith("*"):
                fields = line.split()
                if len(fields) >= 10 and fields[2] == station and fields[1] in ("STAX", "STAY", "STAZ"):
                    if fields[6] != "m" or fields[1] in axes:
                        raise ValueError("Ambiguous station solution or non-metre SINEX coordinates")
                    axes[fields[1]] = float(fields[8])
                    uncertainties[fields[1]] = float(fields[9])
                    epochs.add(fields[5])
    if len(axes) != 3 or len(epochs) != 1:
        raise ValueError("Missing consistent station XYZ in SOLUTION/ESTIMATE")
    return {"station": station, "epoch": epochs.pop(),
            "xyz_m": [axes[k] for k in ("STAX", "STAY", "STAZ")],
            "formal_sigma_m": [uncertainties[k] for k in ("STAX", "STAY", "STAZ")]}


def enu(vector, xyz):
    # Bowring geodetic latitude, independent from the iterative C++ conversion.
    a, f = 6378137.0, 1/298.257223563
    b, e2 = a*(1-f), f*(2-f)
    p = math.hypot(xyz[0], xyz[1])
    theta = math.atan2(xyz[2]*a, p*b)
    lat = math.atan2(xyz[2]+(a*a-b*b)/(b*b)*b*math.sin(theta)**3,
                     p-e2*a*math.cos(theta)**3)
    lon = math.atan2(xyz[1], xyz[0])
    sl, cl, sp, cp = math.sin(lon), math.cos(lon), math.sin(lat), math.cos(lat)
    x, y, z = vector
    return [-sl*x+cl*y, -sp*cl*x-sp*sl*y+cp*z, cp*cl*x+cp*sl*y+sp*z]


def position_comparison(path, station):
    xyz = station["xyz_m"]
    differences = []
    with open(path, encoding="utf-8-sig", newline="") as source:
        for r in csv.DictReader(source):
            position = [float(r[k]) for k in ("x", "y", "z")]
            if all(map(math.isfinite, position)):
                differences.append(enu([a-b for a, b in zip(position, xyz)], xyz))
    return {"epochs": len(differences), "enu_m": {
        key: stats(v[i] for v in differences) for i, key in enumerate(("E", "N", "U"))},
        "three_dimensional_m": stats(norm(v) for v in differences)}


def sp3_velocity(records, sat, time, step):
    # Ninth-point centered derivative. Drop four epochs at each end; no extrapolation.
    pairs = [(4/5, 1), (-1/5, 2), (4/105, 3), (-1/280, 4)]
    if any((sat, time+s*step*k) not in records for _, k in pairs for s in (-1, 1)):
        return None
    return [sum(w*(records[(sat, time+step*k)][axis]-records[(sat, time-step*k)][axis])
                for w, k in pairs)/step for axis in range(3)]


def checksum(path):
    h = hashlib.sha256()
    with open(path, "rb") as source:
        for chunk in iter(lambda: source.read(1024*1024), b""):
            h.update(chunk)
    return {"filename": Path(path).name, "sha256": h.hexdigest()}


def compare(args):
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    precise, epochs, sp3_header = load_sp3(args.sp3)
    clocks, clk_header = load_clk(args.clk)
    station = load_sinex(args.sinex, args.station)
    step = epochs[1]-epochs[0]
    if any(abs(b-a-step) > 1e-6 for a, b in zip(epochs, epochs[1:])):
        raise ValueError("SP3 epochs must have a constant interval")
    queries = sorted(precise, key=lambda k: (k[1], k[0]))
    text = "".join(f"{int(time//604800)} {time%604800:.6f} {sat}\n" for sat, time in queries)
    command = [str(Path(args.state_dump).resolve()), str(Path(args.nav).resolve())]
    run = subprocess.run(command, input=text, capture_output=True, text=True, check=True)
    (output/"broadcast_states.csv").write_text(run.stdout, encoding="utf-8")
    rows, unavailable = [], 0
    for state in csv.DictReader(io.StringIO(run.stdout)):
        if state["status"] != "ok":
            unavailable += 1
            continue
        sat = state["sat"]
        time = int(state["week"])*604800+float(state["tow"])
        p = precise[(sat, time)]
        delta = [float(state[k])-p[i] for i, k in enumerate(("x", "y", "z"))]
        row = {"sat": sat, "week": int(state["week"]), "tow": float(state["tow"]),
               "orbit_3d_m": norm(delta), "radial_m": None, "along_m": None, "cross_m": None,
               "velocity_3d_mps": None, "clock_raw_m": None, "clock_relative_m": None,
               "drift_raw_mps": None, "drift_relative_mps": None}
        v = sp3_velocity(precise, sat, time, step)
        if v:
            vi = [v[0]-OMEGA*p[1], v[1]+OMEGA*p[0], v[2]]
            radial, normal = unit(p), unit(cross(p, vi))
            along = cross(normal, radial)
            row.update(radial_m=dot(delta, radial), along_m=dot(delta, along), cross_m=dot(delta, normal))
            row["velocity_3d_mps"] = norm([float(state[k])-v[i] for i, k in enumerate(("vx", "vy", "vz"))])
        if (sat, time) in clocks:
            row["clock_raw_m"] = C*(float(state["clock_polynomial"])-clocks[(sat, time)])
        if (sat, time-30) in clocks and (sat, time+30) in clocks:
            precise_drift = (clocks[(sat, time+30)]-clocks[(sat, time-30)])/60
            row["drift_raw_mps"] = C*(float(state["drift_polynomial"])-precise_drift)
        rows.append(row)
    if not rows:
        raise ValueError("No usable broadcast states matched the precise orbit epochs")
    groups = defaultdict(list)
    for row in rows:
        groups[(row["sat"][0], row["week"], row["tow"])].append(row)
    offsets = defaultdict(list)
    for (sys, _, _), group in groups.items():
        for raw, relative in (("clock_raw_m", "clock_relative_m"), ("drift_raw_mps", "drift_relative_mps")):
            valid = [r for r in group if r[raw] is not None]
            if len(valid) < 2:
                continue
            offset = statistics.mean(r[raw] for r in valid)
            offsets[(sys, raw)].append(offset)
            for r in valid:
                r[relative] = r[raw]-offset
    metrics = ("orbit_3d_m", "radial_m", "along_m", "cross_m", "velocity_3d_mps",
               "clock_raw_m", "clock_relative_m", "drift_raw_mps", "drift_relative_mps")
    summary = {
        "products": {k: checksum(getattr(args, k)) for k in ("sp3", "clk", "sinex", "nav")},
        "sp3_header": sp3_header, "clk_header_notes": [s for s in clk_header if "COMMENT" in s],
        "station": station, "sp3_epoch_count": len(epochs), "sp3_interval_s": step,
        "queried_satellite_epochs": len(queries), "unavailable_broadcast_epochs": unavailable,
        "systems": {}, "worst_satellites_by_orbit_rms": [], "position_vs_sinex": {},
    }
    for sys in "GCER":
        subset = [r for r in rows if r["sat"][0] == sys]
        summary["systems"][sys] = {key: stats(r[key] for r in subset if r[key] is not None) for key in metrics}
        summary["systems"][sys]["epoch_clock_common_offset_m"] = stats(offsets[(sys, "clock_raw_m")])
    for sat in sorted({r["sat"] for r in rows}):
        subset = [r for r in rows if r["sat"] == sat]
        summary["worst_satellites_by_orbit_rms"].append({"sat": sat, **stats(r["orbit_3d_m"] for r in subset)})
    summary["worst_satellites_by_orbit_rms"].sort(key=lambda r: r["rms"], reverse=True)
    summary["worst_satellites_by_orbit_rms"] = summary["worst_satellites_by_orbit_rms"][:12]
    for path in args.positions:
        summary["position_vs_sinex"][Path(path).name] = position_comparison(path, station)
    with (output/"differences.csv").open("w", encoding="utf-8", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output/"summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("sp3", "clk", "sinex", "nav", "state-dump", "output"):
        parser.add_argument("--"+name, required=True)
    parser.add_argument("--station", default="JFNG")
    parser.add_argument("--positions", nargs="*", default=[])
    summary = compare(parser.parse_args())
    print(json.dumps({"station": summary["station"], "systems": summary["systems"],
                      "position_vs_sinex": summary["position_vs_sinex"]}, indent=2))

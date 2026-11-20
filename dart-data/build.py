"""Build the DART window tracker's data from DART's published GTFS timetable.

Usage: python build.py <path-to-google_transit.zip-or-folder> <out-dir>
Writes rail.json (track shapes, stations, run-time profiles) and sched.json
(every rail trip with its stop times and distances along its track shape).
"""
import csv, io, json, math, os, statistics, sys, zipfile, collections

src, out = sys.argv[1], sys.argv[2]
if os.path.isdir(src):
    opener = lambda f: open(os.path.join(src, f), encoding="utf-8-sig")
else:
    z = zipfile.ZipFile(src)
    names = {n.split("/")[-1]: n for n in z.namelist()}
    opener = lambda f: io.TextIOWrapper(z.open(names[f]), encoding="utf-8-sig")
rd = lambda f: csv.DictReader(opener(f))

RAIL_TYPES = ("0", "1", "2", "5")
routes = {r["route_id"]: r for r in rd("routes.txt") if r["route_type"] in RAIL_TYPES}
trips = {t["trip_id"]: t for t in rd("trips.txt") if t["route_id"] in routes}
shape_ids = {t["shape_id"] for t in trips.values()}

pts = collections.defaultdict(list)
for p in rd("shapes.txt"):
    if p["shape_id"] in shape_ids:
        pts[p["shape_id"]].append((int(p["shape_pt_sequence"]), float(p["shape_pt_lat"]), float(p["shape_pt_lon"])))

def tsec(s):
    h, m, x = s.strip().split(":")
    return int(h) * 3600 + int(m) * 60 + int(x)

st = collections.defaultdict(list)
for r in rd("stop_times.txt"):
    if r["trip_id"] in trips and r["arrival_time"].strip():
        st[r["trip_id"]].append((int(r["stop_sequence"]), r["stop_id"], tsec(r["arrival_time"]), tsec(r["departure_time"] or r["arrival_time"])))
stops_needed = {x[1] for v in st.values() for x in v}
stops = {s["stop_id"]: s for s in rd("stops.txt") if s["stop_id"] in stops_needed}

LAT0 = 32.8
KX = 111320 * math.cos(math.radians(LAT0)); KY = 110540
def xy(lat, lon): return ((lon + 96.8) * KX, (lat - LAT0) * KY)

def simplify(P, tol):
    if len(P) < 3: return P
    keep = [False] * len(P); keep[0] = keep[-1] = True
    stack = [(0, len(P) - 1)]
    while stack:
        a, b = stack.pop()
        ax, ay = P[a][2]; bx, by = P[b][2]; dx, dy = bx - ax, by - ay; L = math.hypot(dx, dy) or 1e-9
        best, idx = -1, None
        for i in range(a + 1, b):
            px, py = P[i][2]; d = abs(dy * px - dx * py + bx * ay - by * ax) / L
            if d > best: best, idx = d, i
        if best > tol:
            keep[idx] = True; stack += [(a, idx), (idx, b)]
    return [p for p, k in zip(P, keep) if k]

shapes = {}
for sid, L in pts.items():
    L.sort()
    shapes[sid] = simplify([(la, lo, xy(la, lo)) for _, la, lo in L], 4.0)

def cumlen(S):
    c = [0.0]
    for i in range(1, len(S)):
        c.append(c[-1] + math.hypot(S[i][2][0] - S[i - 1][2][0], S[i][2][1] - S[i - 1][2][1]))
    return c
CUM = {sid: cumlen(S) for sid, S in shapes.items()}

def project(sid, x, y, after=0.0):
    """Distance along shape of the closest point to (x,y), searching forward from `after`."""
    S, C = shapes[sid], CUM[sid]
    best = (1e18, 0.0)
    for i in range(len(S) - 1):
        if C[i + 1] < after - 5: continue
        ax, ay = S[i][2]; bx, by = S[i + 1][2]; dx, dy = bx - ax, by - ay; L2 = dx * dx + dy * dy or 1e-9
        t = max(0, min(1, ((x - ax) * dx + (y - ay) * dy) / L2))
        d = math.hypot(x - ax - t * dx, y - ay - t * dy)
        if d < best[0]: best = (d, C[i] + t * math.sqrt(L2))
    return best[1]

# stop distances along each shape, per stop pattern (cached)
pat_cache = {}
def pattern_dists(sid, stop_ids):
    key = (sid, stop_ids)
    if key not in pat_cache:
        ds, after = [], 0.0
        for s in stop_ids:
            x, y = xy(float(stops[s]["stop_lat"]), float(stops[s]["stop_lon"]))
            d = project(sid, x, y, after)
            if d < after: d = after
            ds.append(round(d)); after = d
        pat_cache[key] = ds
    return pat_cache[key]

shape_trips = collections.defaultdict(list)
for tid, t in trips.items(): shape_trips[t["shape_id"]].append(tid)

# ---- rail.json: shapes, stations, a typical run-time profile per shape ----
out_shapes = []
for sid, S in shapes.items():
    tl = [t for t in shape_trips[sid] if st[t]]
    if not tl: continue
    t0 = trips[tl[0]]; r = routes[t0["route_id"]]
    head = collections.Counter(trips[x]["trip_headsign"] for x in tl).most_common(1)[0][0]
    pattern = collections.Counter(tuple(x[1] for x in sorted(st[t])) for t in tl).most_common(1)[0][0]
    offs = collections.defaultdict(list)
    for t in tl:
        v = sorted(st[t])
        if tuple(x[1] for x in v) != pattern: continue
        base = v[0][3]
        for i, x in enumerate(v): offs[i].append((x[2] - base, x[3] - base))
    ds = pattern_dists(sid, pattern)
    prof = [[ds[i], round(statistics.median(o[0] for o in offs[i])), round(statistics.median(o[1] for o in offs[i])), s] for i, s in enumerate(pattern)]
    out_shapes.append({"id": sid, "route": r["route_short_name"], "route_id": r["route_id"], "color": "#" + (r.get("route_color") or "888888"),
                       "type": int(r["route_type"]), "head": head, "dir": t0.get("direction_id", ""), "ntrips": len(tl),
                       "pts": [[round(la, 5), round(lo, 5)] for la, lo, _ in S], "prof": prof})

feed_version = ""
try:
    fi = next(rd("feed_info.txt")); feed_version = fi.get("feed_version", "")
except Exception: pass

rail = {"version": feed_version,
        "routes": {k: [v["route_short_name"], "#" + (v.get("route_color") or "888888"), v["route_long_name"]] for k, v in routes.items()},
        "shapes": out_shapes,
        "stops": {s: [v["stop_name"], round(float(v["stop_lat"]), 5), round(float(v["stop_lon"]), 5)] for s, v in stops.items()}}

# ---- sched.json: every rail trip ----
services = {}
try:
    for c in rd("calendar.txt"):
        services[c["service_id"]] = [[int(c[d]) for d in ("monday", "tuesday", "wednesday", "thursday", "friday", "saturday", "sunday")], int(c["start_date"]), int(c["end_date"])]
except KeyError: pass
exceptions = collections.defaultdict(dict)
try:
    for c in rd("calendar_dates.txt"):
        exceptions[c["date"]][c["service_id"]] = int(c["exception_type"])
except KeyError: pass

used_services = set()
trip_rows = []
for tid, t in trips.items():
    v = sorted(st[tid])
    if len(v) < 2 or t["shape_id"] not in shapes: continue
    ds = pattern_dists(t["shape_id"], tuple(x[1] for x in v))
    start = v[0][3]
    # [shape_id, service_id, start_sec, [dist...], [arr offset...], [dep offset...]]
    trip_rows.append([t["shape_id"], t["service_id"], start, ds, [x[2] - start for x in v], [x[3] - start for x in v]])
    used_services.add(t["service_id"])
trip_rows.sort(key=lambda r: (r[1], r[2]))

sched = {"version": feed_version, "tz": "America/Chicago",
         "services": {k: v for k, v in services.items() if k in used_services},
         "exceptions": {d: {s: e for s, e in m.items() if s in used_services} for d, m in exceptions.items() if any(s in used_services for s in m)},
         "trips": trip_rows}

os.makedirs(out, exist_ok=True)
json.dump(rail, open(os.path.join(out, "rail.json"), "w"), separators=(",", ":"))
json.dump(sched, open(os.path.join(out, "sched.json"), "w"), separators=(",", ":"))
# script versions so the page can load them directly
open(os.path.join(out, "rail.js"), "w").write("window.DART_RAIL=" + json.dumps(rail, separators=(",", ":")) + ";")
open(os.path.join(out, "sched.js"), "w").write("window.DART_SCHED=" + json.dumps(sched, separators=(",", ":")) + ";")
print(f"version {feed_version}: {len(out_shapes)} shapes, {len(trip_rows)} rail trips, {len(sched['services'])} services")

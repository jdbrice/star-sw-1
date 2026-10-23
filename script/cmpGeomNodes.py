#!/usr/bin/env python3
"""Difference two dumpGeomNodes.C outputs: misaligned - ideal, per node.

Groups by detector (STFM = sTGC pentagon, FTUS = FST sensor) and reports every
distinct shift, so a geometry that moved in a few rigid groups shows up as a few
rows rather than 124.
"""
import sys, collections

def read(fn):
    d = {}
    for line in open(fn):
        p = line.split()
        if len(p) != 4: continue
        d[p[0]] = (float(p[1]), float(p[2]), float(p[3]))
    return d

ideal, misal = read(sys.argv[1]), read(sys.argv[2])
only_i = set(ideal) - set(misal)
only_m = set(misal) - set(ideal)
if only_i or only_m:
    print("WARNING: node sets differ -- %d only ideal, %d only misaligned" % (len(only_i), len(only_m)))

groups = collections.defaultdict(list)
for path in sorted(set(ideal) & set(misal)):
    a, b = ideal[path], misal[path]
    d = tuple(round(b[k] - a[k], 4) for k in range(3))
    det = "STFM" if "STFM" in path else ("FTUS" if "FTUS" in path else "other")
    groups[(det, d)].append(path)

print("%-6s %10s %10s %10s %6s   %s" % ("det", "dx(cm)", "dy(cm)", "dz(cm)", "nodes", "example"))
tot = 0
for (det, d), paths in sorted(groups.items(), key=lambda kv: (kv[0][0], -len(kv[1]))):
    tot += len(paths)
    ex = paths[0].split('/')[-1]
    print("%-6s %10.4f %10.4f %10.4f %6d   %s%s" %
          (det, d[0], d[1], d[2], len(paths), ex, "" if len(paths) == 1 else " ..."))
print("\n%d nodes compared, %d distinct shifts" % (tot, len(groups)))
moved = sum(len(p) for (det, d), p in groups.items() if any(abs(v) > 1e-4 for v in d))
print("%d moved, %d did not" % (moved, tot - moved))

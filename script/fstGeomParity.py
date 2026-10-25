#!/usr/bin/env python3
"""Front/back mounting parity of every FST wedge, from a dumpGeomNodes.C dump.

FTUS_<n> with n = 1 + 3*iwedgem + 36*idisk + sensor, so
   idisk = (n-1)//36,  iwedgem = ((n-1)%36)//3,  sensor = (n-1)%3.
A wedge is mounted FRONT (silicon at z < disk nominal) or BACK (z > nominal);
AGML drives that with wholewedge_modulus2, which fstmconfig selects.
"""
import sys, collections
ZDISK = (151.750, 165.248, 178.781)

def read(fn):
    out = {}
    for line in open(fn):
        p = line.split()
        if len(p) != 4 or not p[0].startswith("FTUS_"):
            continue
        n = int(p[0].split("_")[1])
        out[n] = float(p[3])
    return out

for fn in sys.argv[1:]:
    z = read(fn)
    print("=== %s : %d FTUS nodes" % (fn, len(z)))
    for idisk in range(3):
        row = []
        for iw in range(12):
            n = 1 + 3*iw + 36*idisk          # inner sensor of this wedge
            if n not in z: row.append(" ? "); continue
            dz = z[n] - ZDISK[idisk]
            row.append("%+d" % (1 if dz > 0 else -1))
        print("  disk %d (nominal z %8.3f): %s" % (idisk+1, ZDISK[idisk], " ".join(row)))
    # the actual offsets, one per (disk, sensor class)
    off = collections.defaultdict(set)
    for n, zz in z.items():
        idisk = (n-1)//36; sensor = (n-1) % 3
        off[(idisk, sensor)].add(round(abs(zz - ZDISK[idisk]), 4))
    print("  |dz| per (disk, sensor):")
    for k in sorted(off):
        print("    disk %d sensor %d : %s" % (k[0]+1, k[1], sorted(off[k])))

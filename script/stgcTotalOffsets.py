#!/usr/bin/env python3
"""Total sTGC offset per station/quadrant = what StFttDb actually applies to a hit.

    total = AGML placeholder + stgcOnTpc + stationOnStgc + pentOnStation

matching StFttDb::loadGeometryFromDb(). AGML placeholder is x = {0,0,-6.5,+6.5} by
PENTAGON index and y = +5.9 for all (the STGM mother offset in StgmGeo1.xml).
Pentagon order within a station is A,D,C,B, i.e. pent = (4-quad)%4, so row = 4*st+pent.

Usage: stgcTotalOffsets.py <stgcOnTpc.C> <stationOnStgc.C> <pentOnStation.C>
"""
import re, sys

AGML_X = [0.0, 0.0, -6.5, 6.5]     # by pentagon (A, D, C, B)
AGML_Y = 5.9
QUAD2PENT = {"A": 0, "B": 3, "C": 2, "D": 1}

def rows(path):
    """[(t0,t1,t2), ...] in file order."""
    txt = open(path, encoding="utf-8", errors="replace").read()
    out = []
    for m in re.finditer(r'row\.t0\s*=\s*([-0-9.eE]+);\s*\n\s*row\.t1\s*=\s*([-0-9.eE]+);\s*\n\s*row\.t2\s*=\s*([-0-9.eE]+);', txt):
        out.append(tuple(float(g) for g in m.groups()))
    return out

def totals(f_tpc, f_station, f_pent):
    g = rows(f_tpc)[0] if rows(f_tpc) else (0.0, 0.0, 0.0)
    st = rows(f_station)
    pe = rows(f_pent)
    res = {}
    for s in range(4):
        for q in "ABCD":
            p = QUAD2PENT[q]
            r = 4*s + p
            sx, sy, sz = st[s] if s < len(st) else (0.0, 0.0, 0.0)
            px, py, pz = pe[r] if r < len(pe) else (0.0, 0.0, 0.0)
            res[(s, q)] = (AGML_X[p] + g[0] + sx + px,
                           AGML_Y    + g[1] + sy + py)
    return res

if __name__ == "__main__":
    t = totals(*sys.argv[1:4])
    for s in range(4):
        for q in "ABCD":
            print("  station %d quad %s   x = %+8.4f   y = %+8.4f  cm" % (s+1, q, t[(s,q)][0], t[(s,q)][1]))

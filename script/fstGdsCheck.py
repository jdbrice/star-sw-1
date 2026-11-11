#!/usr/bin/env python3
"""Read the FST sensor mask files and check our geometry constants against them.

WHY THIS EXISTS
    The outer-overlap measurement fits two parameters -- a fractional strip
    PITCH error and a RADIAL OFFSET of the point the strip fan converges on --
    and both are properties of the silicon, so both can be checked against the
    thing the silicon was cut from. They do not survive the check: the mask
    gives the pitch as exactly 30/128 deg and puts every strip line within
    0.1 um of the beam axis. Whatever widens the measured overlap is in how the
    wedges are MOUNTED, not in the sensors.

    Inputs (akio's copies, read-only):
        ~/www/fcs/survey/STAR-FST-INNER-190301.gds
        ~/www/fcs/survey/STAR-FST-OUTER-190313.gds
    and the design document Speci-V.005.1.pdf alongside them.

GDSII, the little of it that is needed
    A stream of records: 2-byte big-endian length (header included), 1-byte
    record type, 1-byte data type, then the payload.
      0x0305 UNITS    two 8-byte GDS reals: user units per db unit, METRES per
                      db unit. Here 1 db unit = 1 nm.
      0x0606 STRNAME  the structure being defined
      0x0D02 LAYER
      0x1003 XY       int32 pairs, in db units
      0x0A00 SREF / 0x0B00 AREF  a placement; 0x1206 SNAME names what is
                      placed, 0x1C05 ANGLE is its rotation in degrees
      0x1100 ENDEL    ends the element
    A GDS real is sign bit, 7-bit excess-64 exponent base 16, 56-bit fraction.

Usage:
    python3 -I script/fstGdsCheck.py [inner.gds outer.gds]
"""

import math
import os
import struct
import sys


def gds_real(b):
    sign = -1.0 if b[0] & 0x80 else 1.0
    exp = (b[0] & 0x7F) - 64
    frac = int.from_bytes(b[1:8], 'big') / float(1 << 56)
    return sign * frac * (16.0 ** exp)


def read_gds(path):
    """Return (placements, polygons, cm_per_db_unit)."""
    data = open(path, 'rb').read()
    i = 0
    dbu = None
    sname = None
    layer = None
    placements = []
    polygons = []
    cur = None
    while i + 4 <= len(data):
        length = struct.unpack('>H', data[i:i + 2])[0]
        if length < 4:
            break
        rtype, dtype = data[i + 2], data[i + 3]
        body = data[i + 4:i + length]
        if rtype == 0x03 and dtype == 0x05:
            dbu = gds_real(body[8:16])
        elif rtype == 0x06:
            sname = body.rstrip(b'\x00').decode('ascii', 'replace')
        elif rtype == 0x0D:
            layer = struct.unpack('>h', body[0:2])[0]
        elif rtype in (0x0A, 0x0B):
            cur = {'in': sname, 'ref': None, 'angle': 0.0, 'xy': []}
        elif rtype == 0x12 and cur is not None:
            cur['ref'] = body.rstrip(b'\x00').decode('ascii', 'replace')
        elif rtype == 0x1C and cur is not None:
            cur['angle'] = gds_real(body[0:8])
        elif rtype == 0x10 and dtype == 0x03:
            pts = [struct.unpack('>ii', body[8 * k:8 * k + 8])
                   for k in range(len(body) // 8)]
            if cur is not None:
                cur['xy'] = pts
            else:
                polygons.append((sname, layer, pts))
        elif rtype == 0x11 and cur is not None:
            placements.append(cur)
            cur = None
        i += length
    return placements, polygons, dbu * 100.0


def check(path, label, r_cut, expect_half_fan_deg):
    placements, polygons, sc = read_gds(path)

    # The sensor strips are HPK_strip placed once per (phi strip, r strip).
    # Test structures ("baby" cells) sit at small radius and are cut away.
    strips = []
    for pl in placements:
        if pl['ref'] != 'HPK_strip' or not pl['xy'] or pl['angle'] == 0.0:
            continue
        x, y = pl['xy'][0][0] * sc, pl['xy'][0][1] * sc
        if math.hypot(x, y) < r_cut:
            continue
        strips.append((pl['angle'], x, y))

    print('=== %s  (%s)' % (label, os.path.basename(path)))
    print('  %d sensor strip placements' % len(strips))

    # Distance from the beam axis to each strip LINE. A strip is a bar through
    # (x, y) in the direction of its own rotation angle, so the perpendicular
    # distance from the origin is |x sin a - y cos a|. If the fan converges on
    # the axis this is zero for every strip.
    dists = [(x * math.sin(math.radians(a)) - y * math.cos(math.radians(a))) * 1e4
             for a, x, y in strips]                       # microns
    print('  fan vertex: every strip line passes within %.4f um of the beam axis'
          % max(abs(d) for d in dists))
    print('              (mean %+.5f um, min %+.5f, max %+.5f)'
          % (sum(dists) / len(dists), min(dists), max(dists)))

    angles = sorted(set(round(a, 5) for a, _, _ in strips))
    steps = [angles[k + 1] - angles[k] for k in range(len(angles) - 1)]
    print('  pitch: %d distinct angles, step %.6f to %.6f deg  (30/128 = %.6f)'
          % (len(angles), min(steps), max(steps), 30.0 / 128))
    # Only half the fan is drawn; the other half is mirrored about the sensor
    # bisector at 90 deg. First strip CENTRE to last strip centre is therefore
    # 2*(90 - angles[0]), and the sensor spans that plus one pitch of edges.
    full = 2 * (90.0 - angles[0]) + 30.0 / 128
    print('  fan: drawn %.5f..%.5f deg -> full sensor %.6f deg  (expected %.3f)'
          % (angles[0], angles[-1], full, expect_half_fan_deg))

    radii = sorted(set(round(math.hypot(x, y), 4) for _, x, y in strips))
    print('  strip-start radii: %s cm' % radii)

    # The die outline: the substrate polygon with the largest radial span. Its
    # straight side edge is NOT radial -- it misses the axis by a fixed amount,
    # which makes the die a trapezoid slightly wider than its sector. That is
    # dicing and guard-ring margin, outside the active area.
    best = None
    for sname, lay, pts in polygons:
        if not (3 <= len(pts) <= 60):
            continue
        rs = [math.hypot(x * sc, y * sc) for x, y in pts]
        if max(rs) < 3.0:
            continue
        if best is None or (max(rs) - min(rs)) > best[0]:
            best = (max(rs) - min(rs), sname, lay, pts)
    if best:
        _, sname, lay, pts = best
        # The SIDE edge runs between the two vertices furthest off the
        # bisector -- one on the inner arc, one on the outer. Taking min-r and
        # max-r instead picks two points on different edges and is meaningless.
        vs = [(x * sc, y * sc) for x, y in pts]
        rr_all = [math.hypot(x, y) for x, y in vs]
        rmid = 0.5 * (min(rr_all) + max(rr_all))
        off = lambda v: abs(math.degrees(math.atan2(v[0], v[1])))
        lo = [v for v in vs if math.hypot(*v) < rmid]
        hi = [v for v in vs if math.hypot(*v) >= rmid]
        # one endpoint on the inner arc and one on the outer, both the furthest
        # off the bisector: taking the two largest overall picks both ends of
        # the same arc, because the outline is drawn as two concentric contours.
        (x1, y1) = max(lo, key=off)
        (x2, y2) = max(hi, key=off)
        dx, dy = x2 - x1, y2 - y1
        perp = abs(x1 * dy - y1 * dx) / math.hypot(dx, dy) * 10.0     # mm
        rr = [math.hypot(x, y) for x, y in vs]
        print('  die outline (%s layer %s): r %.4f to %.4f cm' %
              (sname, lay, min(rr), max(rr)))
        print('    side edge from r = %.4f to %.4f cm misses the axis by '
              '%.4f mm, so the die is a trapezoid a little wider than its '
              'sector -- dicing margin, outside the active area'
              % (math.hypot(x1, y1), math.hypot(x2, y2), perp))
    print()


def main(argv):
    home = os.path.expanduser('~')
    inner = argv[1] if len(argv) > 2 else \
        home + '/www/fcs/survey/STAR-FST-INNER-190301.gds'
    outer = argv[2] if len(argv) > 2 else \
        home + '/www/fcs/survey/STAR-FST-OUTER-190313.gds'
    check(inner, 'INNER', 4.0, 30.0)
    check(outer, 'OUTER', 10.0, 15.0)
    print('Our constants, for comparison:')
    print('  StFstConsts.h  kFstStripPitchPhi = 30 deg / 128      = %.6f deg'
          % (30.0 / 128))
    print('                 kFstStripPitchR   = 2.875 cm')
    print('                 kFstrStart = 5.000 7.875 10.750 13.625 '
          '16.500 19.375 22.250 25.125 cm')
    print('  FstmGeoSM.xml  innersilicon rmin=5.0  rmax=16.5  phi -15..+15')
    print('                 outersilicon rmin=16.5 rmax=28.0 '
          'phi {-15.5,-0.5} and {0.5,15.5}')


main(sys.argv)

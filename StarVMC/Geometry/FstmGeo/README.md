# FSTM geometry model variants

This directory contains four related geometry descriptions for the Forward Silicon Tracker (FST). They all describe the same detector family, but they are aimed at different use cases:

- Full, detailed detector description for realistic material and geometry studies
- Simplified detector description for faster simulation and idealized studies
- Misalignment-capable variants for position and mirror-configuration studies

## Files in this directory

### FstmGeo.xml
The original full-geometry model. This is the baseline detailed description of the FST support structure and silicon modules. It includes the full set of mechanical support pieces, ribs, heat sinks, cooling structures, and hybrid/T-board elements, and it retains the full disk/wedge description (`FSTD`, `FSTW`, various rib and heat-sink components, cooling-tube holders, and connectors).

This file is the reference for the detailed detector layout and the geometry from which the simplified versions were derived.

### FstmGeoM.xml
A later full-geometry variant, intended as the full model in the same family but with additional bookkeeping for misalignment and mirror/wedge configuration. It keeps the detailed detector structure, but it adds variables such as `idisk`, `iwedgem`, `halfid`, and wedge/silicon offset parameters for placement control. It also adds the `FSTG.fstmconfig` configuration variable and the 2026 steering support for wedge mirror selection.

Key features:
- Full geometry retained
- Detailed support/rib/cooling structures retained
- Added placement and misalignment controls
- Added mirror-config steering for wedge orientation

The comments in the file also note important later updates, including the 2021 ring relocation and the 2026 mirror-choice steering.

### FstmGeoS.xml
A simplified version of the original full geometry, copied from `FstmGeo.xml` and reduced for speed and cleaner idealized studies. It is still a disk-plus-wedge model, but many small and repetitive components are removed or collapsed.

The simplified model keeps the core detector structure and cooling channels, while eliminating most of the extra rib, heat-sink, and support-holder elements. It replaces the detailed cooling network with a simpler set of straight cooling tubes and combines stainless steel + NOVEC into a single material (`SteelNOVECMix`).

Important simplifications:
- Removed or collapsed outer and inner rib structures
- Removed extra heat-sink pieces and cooling-tube holder elements
- Replaced detailed cooling geometry with a simpler straight-tube model
- Merged material components into a single combined material
- Keeps the detailed disk container (`FSTD`) and wedge assembly structure

### FstmGeoSM.xml
The simplified misaligned geometry variant. This is the most reduced and configuration-aware version of the family. It is effectively the simplified model with additional support for misalignment and wedge-mirror placement logic.

Like `FstmGeoS.xml`, it removes most of the detailed mechanical subcomponents. However, unlike the simplified ideal version, it includes explicit misalignment bookkeeping and later installation/mirror corrections. It also updates the support-ring and wedge-position treatment for the actual installation geometry.

Key additions relative to `FstmGeoS.xml`:
- `codeversion` set to 1.3
- Added `fstmconfig`-style steering control for wedge mirror configuration
- Adjusted wedge placement to reflect real installation geometry
- Retains simplified cooling and merged material model
- Designed for misalignment and configuration studies rather than full-mechanics realism

## High-level comparison

### Full vs. simplified
The crucial split is between:

- `FstmGeo.xml` and `FstmGeoM.xml`: full geometry, detailed mechanical layout
- `FstmGeoS.xml` and `FstmGeoSM.xml`: simplified geometry, reduced volume count and merged materials

### Idealized vs. misaligned
Within each family:

- `FstmGeo.xml` / `FstmGeoS.xml`: baseline full or simplified geometry
- `FstmGeoM.xml` / `FstmGeoSM.xml`: versions with explicit misalignment and mirror configuration support

### Common design choices across simplified models
The simplified files share these major design changes, documented directly in the comments:

- Mass from PEEK inner and outer ribs folded into respective PEEK bases
- Cooling tubes replaced with four straight tubes
- Stainless steel + NOVEC combined into one material: `SteelNOVECMix`
- Cooling-tube holder removed
- Additional cooling-tube connector pieces added in the simplified version
- Heat sinks compressed into a smaller number of volumes inside inner/outer PEEK bases
- Significant commented-out legacy code removed

## Practical interpretation

Use the full geometry files when you need:
- realistic detector material budget
- detailed mechanical support and support-system studies
- precise geometry matching to drawings and installation constraints

Use the simplified geometry files when you need:
- faster simulation setup
- idealized detector studies
- misalignment or mirror-configuration scans
- reduced detector complexity for performance-sensitive workflows

## Summary

In short:

- `FstmGeo.xml` = original detailed full geometry
- `FstmGeoM.xml` = later detailed full geometry with additional placement/misalignment/mirror controls
- `FstmGeoS.xml` = simplified copy of the detailed model with merged materials and reduced component count
- `FstmGeoSM.xml` = simplified misaligned geometry, most reduced and most configuration-aware of the family

The overall evolution is from a detailed mechanical model to a simplified, faster-to-simulate geometry, while also adding explicit control over misalignment and wedge orientation in the later variants.

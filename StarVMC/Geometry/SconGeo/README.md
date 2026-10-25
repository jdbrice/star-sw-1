# SconGeo model comparison

This directory contains two versions of the STAR support-structure geometry for the silicon tracker environment:

- `SconGeo.xml`: the established support-structure model used across much of the STAR lifetime
- `SconGeo15.xml`: a refined, updated geometry that reflects a better understanding of the materials and their distribution

## 1. General purpose

Both files model the support cone and beampipe-support system installed around the silicon tracker in the cave. The geometry includes carbon-composite support rods, support cones, kinematic mounts, roller assemblies, and associated support hardware used to hold the detector system in place.

The two models are complementary rather than mutually exclusive: `SconGeo.xml` is the long-standing support geometry that has been maintained and used extensively, while `SconGeo15.xml` is the updated version designed to better represent the actual material layout and support structure.

## 2. Relationship between the two files

`SconGeo.xml` was originally derived from the historical SVT support model and was later split into a dedicated support geometry file. It remains an actively supported AGML geometry.

`SconGeo15.xml` was produced as a refreshed version of that same family of support models. The file comments indicate that it was written as a cleaner AGML implementation and reflects an improved understanding of materials and the geometry distribution. The newer version is not a replacement of the older one for all use cases; it is a refinement aimed at more accurate description of the support structure.

## 3. Structural similarities

Both models describe the same general detector support system and include the same broad physical features:

- support cone around the beam pipe / SVT region
- carbon-composite support rods
- Nomex and carbon-fiber layering in support elements
- kinematic mount / beam-pipe support ring
- roller supports and support shafts
- mechanical support hardware around the beam pipe

The geometry is organized around a support cone mother volume and a beampipe support region, with nested volumes for the cone shell, support ring, and support hardware.

## 4. Structural differences

### SconGeo.xml
`SconGeo.xml` contains a broad set of historical support blocks and support volumes, including:

- `SCON` / `SCOM` support-cone mother structure
- `SGRA` graphite/epoxy cone
- `SCMY` mylar wrapping
- `SBSP` beampipe support region
- `SAKM` aluminum kinematic mount
- `SCKM` cutout in the kinematic mount
- `SBMM`, `SBMI`, `SBMO` mounting blocks
- `SMRD` aluminum support rod
- `SBRL`, `SBRX` ceramic roller and shaft
- `SBSR`, `SBCR` G10 ring and cutout structure
- `SROD`, `SRON`, `SROI`, `SROH` carbon support rod with Nomex filling and hollow interior

This model is parameter-rich and reflects the established support design used for many years. It includes several support versions and multiple structural fills that encode the geometry evolution of the original support model.

### SconGeo15.xml
`SconGeo15.xml` has the same overall physical intent but a more modernized AGML layout. It includes:

- `SCOM` support-cone mother volume
- `SCON` support cone
- `SNMX` Nomex structure
- `SROD`, `SRON`, `SROI`, `SROH` carbon composite support rod and hollow interior
- `SBSP` beampipe support region
- `SAKM` aluminum support mount
- `SAKH` cutout hole in the kinematic mount
- `SPOK`, `SPOA`, `SPOB`, `SPOC`, `SPOD` spoke assembly supports
- `SASH` aluminum EM shield
- `SDSA`, `SDSK` G10 support ring and disk
- `BOLT` bolts in support hardware
- `SBRL`, `SBRX` roller support and shaft

This version is more explicitly structured around support spokes, shield layers, and G10 support disks, and it reflects a more precise material arrangement in the support system.

## 5. Material and layout updates in SconGeo15.xml

The most important difference is not a change of purpose, but an improved material and distribution model. The newer file explicitly introduces better descriptions of:

- Nomex thickness and support-cone layering
- carbon-shell and Nomex structure around support rods
- additional G10 support disk and ring components
- aluminum screen / shielding layer around the beampipe support
- spoke and bolt hardware associated with the support structure

In other words, the updated geometry emphasizes a more realistic support layout and material distribution rather than a change in detector concept.

## 6. Parameter philosophy

`SconGeo.xml` is a long-lived, parameterized support geometry with a broad set of historical variables and multiple support versions.

`SconGeo15.xml` keeps the same basic parameter style but reorganizes the support model to better reflect the actual support materials and their placement. It is not a radical redesign; it is a refinement of the same support structure with a more current understanding of the detector hardware.

## 7. Summary

The essential difference is this:

- `SconGeo.xml` is the established, maintained support geometry used for much of the STAR lifetime.
- `SconGeo15.xml` is a refined AGML implementation of the same support family, with a better material and geometric description of the support cone, beam-pipe support, and surrounding mechanical hardware.

The updated file is therefore best understood as an improved and more detailed support geometry representation, rather than a different detector concept.

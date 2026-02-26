# AstroSeis

MATLAB codes for computing 3-D seismic wavefields in small planetary bodies using the Boundary Element Method (BEM).

## Quick Start

### Homogeneous model

```matlab
AstroSeis BEM_para
```

### Solid body with a liquid core

```matlab
AstroSeis_liquidcore BEM_para_lc
```

### Plotting results

```matlab
AstroSeis_plot BEM_para
AstroSeis_liquidcore_plot BEM_para_lc
```

## Parameter File Format

The parameter file is a plain-text file that controls the simulation. Below is an annotated example (`BEM_para`):

```
phobos_model_20km0.mat       # mesh file name (.mat)
20000 200 1                  # R (mean radius, m), nmesh (~4*nmesh patches), topo_fold (topography scaling)
my_mesh.mat                  # output mesh file name
demo_out_put.mat             # output wavefield file name
6000 3000 3000 200           # vp (m/s), vs (m/s), rho (kg/m^3), Q (attenuation)
500 0.1 0.3                  # nt (time samples), dt (s), f0 (center frequency, Hz)
single 0                     # source type (single or moment), source scale (10^scale)
1 1 1                        # fx, fy, fz (force components, used if source=single)
1 0 0 1 0 1                  # Mxx,Mxy,Mxz,Myy,Myz,Mzz (used if source=moment)
10e3 0 180                   # source depth (m), latitude (deg), longitude (deg)
```

If the mesh file on line 1 already exists, it is loaded directly. Otherwise, the code generates a new mesh using the parameters on line 2 (reference radius, mesh density, and topography scaling).

## Mesh Generation

AstroSeis builds its mesh in two steps:

1. **Reference sphere**: A triangulated sphere of radius `R` is created using `ParticleSampleSphere` and subdivided once.
2. **Topography**: A height value is added to each node of the sphere. The new radial position of each node becomes `r = R + h`, where `h` is the topography height at that node's (colatitude, longitude).

The built-in Phobos example computes the height from spherical harmonic coefficients (Willner et al., PSS 2014), but any source of topography heights can be used.

### To create a mesh for a new body

You need:
- The **mean radius** of the body
- **Topography data** -- height relative to the mean radius as a function of surface position. This can be in any format (spherical harmonic coefficients, a lat/lon height grid, radial height at each point, etc.)

To generate the mesh separately (without running the full simulation):

```matlab
Mesh_generater BEM_para
```

## Examples

Three pre-computed examples are provided in the `examples/` folder:

- **Phobos** -- homogeneous model with topography from spherical harmonics
- **Dimorphos** -- homogeneous model with shape from a triangulated surface model
- **shifted_liquid_core** -- two-layer model with a liquid core

## Contact

Yuan Tian -- ytian159@gmail.com

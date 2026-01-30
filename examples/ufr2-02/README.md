# UFR2-02: Flow Around a Square Cylinder

LES simulation of flow around a square cylinder at Re = 22,000.

## Acknowledgments

The mesh generation script (`generateMeshDict.pl`) and case setup are based on the original work by **Niklas Nordin** for the ERCOFTAC Classic Collection Database.

## Test Case Description

This is ERCOFTAC UFR2-02 / Case 043 from the Classic Collection Database:
- **Geometry**: Square cylinder (D = 0.04 m) in cross-flow
- **Reynolds number**: Re = 22,000 (based on cylinder width and inlet velocity)
- **Inlet velocity**: U = 0.54 m/s
- **Kinematic viscosity**: nu = 9.82e-7 m²/s
- **Domain**: 1.36 m (L) × 0.56 m (W) × 0.392 m (H)

## Reference

Lyn, D.A., Einav, S., Rodi, W., Park, J.-H. (1995)
*"A laser-Doppler velocimetry study of ensemble-averaged characteristics of the turbulent near wake of a square cylinder"*
J. Fluid Mech. 304, 285-319

## Running the Case

```bash
# Build the spaceTimeWindow library first (from repository root)
./Allwmake

# Then run the case
cd examples/ufr2-02
./Allrun
```

## Mesh Generation

The mesh is generated using a Perl script that creates a blockMeshDict:

```bash
./generateMeshDict.pl 0.04 0.392 0.4 1.36 0.56 0.001 > constant/polyMesh/blockMeshDict
```

Parameters:
- D = 0.04 (cylinder diameter)
- H1 = 0.392 (domain height)
- L1 = 0.4 (inlet to cylinder center)
- L2 = 1.36 (total domain length)
- W = 0.56 (domain width)
- delta = 0.001 (cell size near cylinder)

This generates approximately 400,000 cells with refinement near the cylinder walls.

## Using spaceTimeWindow Extraction

To use this case with the spaceTimeWindow library:

1. Uncomment the `extractSubset` function object in `system/controlDict`
2. Run the simulation with extraction enabled
3. Initialize the subset case:
   ```bash
   cd ../ufr2-02-subset
   spaceTimeWindowInitCase -sourceCase ../ufr2-02
   ```
4. Run the reconstruction:
   ```bash
   pimpleFoam
   ```

## Turbulence Model

- LES with Smagorinsky subgrid-scale model
- Cube root volume delta calculation

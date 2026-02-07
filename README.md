# spaceTimeWindow Library

A library for space-time window extraction and reconstruction of LES simulations, developed for OpenFOAM.

**Tested with:** OpenFOAM v2512 (openfoam.com)

> **Trademark Notice:** This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM and OpenCFD trade marks. OPENFOAM is a registered trade mark of OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com.

## Overview

This library provides components for extracting a spatial subset from a full LES simulation and reconstructing the flow within that subset using pre-computed boundary conditions:

1. **spaceTimeWindowExtract** - Function object to extract boundary data during the original simulation
2. **spaceTimeWindowInitCase** - Utility to initialize a reconstruction case from extracted data
3. **pimpleFoamPrecise** - Modified pimpleFoam that registers intermediate velocity HbyA for highest fidelity reconstruction
4. **spaceTimeWindowCoupledPressure** - Pressure BC that blends Dirichlet/Neumann based on flux direction
5. **spaceTimeWindow** / **spaceTimeWindowInletOutlet** - Boundary conditions to apply extracted data during reconstruction

## Key Concepts

### Boundary Condition Approaches

The library provides four approaches for applying boundary conditions on the extraction boundary (`oldInternalFaces`):

#### 1. Precise BC (`-preciseBC`) - **HIGHEST FIDELITY**

Uses the intermediate velocity **u\*** = H(U)/A (called `HbyA` in OpenFOAM) instead of the final velocity U as the boundary condition. Based on Wu, Zaki, Meneveau (Phys. Rev. Fluids 5, 064607, 2020), this is the theoretically optimal BC for subdomain pressure recovery.

**Why u\* is the optimal boundary condition:**
- The pressure Poisson equation is: div(1/A grad(p)) = div(u\*)
- u\* is the *source term* of the pressure equation
- Prescribing u\* exactly on the boundary makes the subdomain pressure equation identical to the full-domain one
- Face values use OpenFOAM's exact distance-weighted interpolation formula — no interpolation error
- Wu et al. (2020) demonstrated machine-precision accuracy in DNS, where the very high resolution minimizes the discretization stencil mismatch at boundary cells (Jasak 1996; Ferziger & Peric 2020), and spectral/direct (non-iterative) pressure solvers yield exact solutions
- On practical LES/RANS meshes with iterative solvers (e.g. GAMG), the coarser grid amplifies this mismatch, though it remains negligible compared to the inherent modelling error of turbulence closures. `-preciseBC` remains the most accurate option (~20% lower pressure errors than standard coupled pressure, see Accuracy Comparison below)

**Requirements:**
- Extraction must use `pimpleFoamPrecise` (registers HbyA in objectRegistry)
- Extraction configuration: `fields (HbyA p)` with `extractGradients true` and `gradientFields (p)`

**Extraction configuration:**
```cpp
extractSubset
{
    type            spaceTimeWindowExtract;
    libs            (spaceTimeWindow);
    box             ((0.05 -0.15 0.01) (0.70 0.15 0.38));
    outputDir       "../subset";

    fields          (HbyA p);       // HbyA from pimpleFoamPrecise + p
    initialFields   (U p nut);      // U (final velocity) for solver initialization

    extractGradients    true;
    gradientFields      (p);

    writeFormat     deltaVarint;
    writeControl    timeStep;
    writeInterval   1;
}
```

**Reconstruction boundary conditions:**
```cpp
// In 0/U (reads HbyA data, applied to U field)
oldInternalFaces
{
    type            spaceTimeWindow;
    dataDir         "constant/boundaryData";
    fieldTableName  HbyA;       // Read HbyA data files
    fixesValue      true;
    allowTimeInterpolation  true;
    timeInterpolationScheme cubic;
    value           uniform (0 0 0);
}

// In 0/p (same as pressure-coupled)
oldInternalFaces
{
    type            spaceTimeWindowCoupledPressure;
    dataDir         "constant/boundaryData";
    blendingMode    fluxMagnitude;
    dirichletWeight 0.8;
    neumannWeight   0.8;
    transitionWidth 0.1;
    allowTimeInterpolation  true;
    timeInterpolationScheme cubic;
    value           uniform 0;
}
```

```bash
# Extraction with pimpleFoamPrecise
cd source-case
pimpleFoamPrecise

# Initialization with -preciseBC
cd ../subset-case
spaceTimeWindowInitCase -sourceCase ../source-case -preciseBC

# Reconstruction with standard pimpleFoam
pimpleFoam
```

#### Accuracy Comparison

Quantitative comparison on the UFR2-02 square cylinder case (110,720 cells full domain, 47,952 cells subset, 88 timesteps with dt=1e-4, no time interpolation, binary uncompressed data):

| Method | U L_inf | U L_2 | p L_inf | p L_2 |
|--------|---------|-------|---------|-------|
| `-preciseBC` (u* + gradp) | 3.71e-3 | 3.07e-4 | 2.19e-3 | 3.22e-4 |
| Coupled pressure (U + gradp) | 3.65e-3 | 3.10e-4 | 2.65e-3 | 3.50e-4 |

Both methods start from identical binary initial conditions (zero error at t_0). Key findings:

- **Velocity errors are nearly identical** — the choice of u* vs U for the velocity BC has negligible effect on velocity reconstruction.
- **Pressure errors are ~20% lower with `-preciseBC`**: L_inf = 2.19e-3 vs 2.65e-3.
- **Errors are concentrated near the boundary.** Only 0.25% of cells exceed 90% of the maximum error. The top-10 worst cells are all at x=0.055 (one cell layer from the extraction box face at x=0.05).

Error decay with distance from the oldInternalFaces boundary (`-preciseBC`):

| Distance from boundary | Cells | U mean err | U max err | p mean err | p max err |
|------------------------|-------|------------|-----------|------------|-----------|
| 0.005–0.01 (1 cell layer) | 1,296 | 1.2e-3 | 3.7e-3 | 1.0e-3 | 2.2e-3 |
| 0.01–0.02 | 5,040 | 2.5e-4 | 9.9e-4 | 4.0e-4 | 1.5e-3 |
| 0.02–0.05 | 12,384 | 1.1e-4 | 6.5e-4 | 2.4e-4 | 8.3e-4 |
| 0.05–0.10 | 15,156 | 4.5e-5 | 2.9e-4 | 1.5e-4 | 7.5e-4 |
| 0.10–0.20 | 14,076 | 1.2e-5 | 9.2e-5 | 5.3e-5 | 5.5e-4 |

In this comparison, no temporal interpolation is used (boundary data at every solver timestep) and face values use OpenFOAM's exact distance-weighted interpolation formula. Despite this, a reconstruction error remains. It arises from the non-linear discretization difference between full-domain internal faces and subdomain Dirichlet boundary faces (Jasak 1996; Moukalled et al. 2016): convective flux limiters (TVD/upwind), the multigrid pressure solver (GAMG), and gradient reconstruction all behave differently at boundary faces (Ferziger & Peric 2020), introducing small differences that compound over timesteps. Wu et al. (2020) report machine-precision accuracy in DNS, where the very high spatial resolution minimizes this stencil mismatch and spectral/direct (non-iterative) pressure solvers yield exact solutions. Note that LES and RANS are inherently approximate---turbulence closure models introduce modelling errors that dwarf the reconstruction error.

#### 2. Pressure-Coupled BC

Uses `spaceTimeWindowCoupledPressure` BC for pressure and `spaceTimeWindow` for velocity. This is the **recommended approach** because it properly handles the elliptic nature of the pressure equation.

**Why pressure coupling is essential:**

Unlike velocity (which is advected and hyperbolic), pressure is governed by an **elliptic equation**. Cutting the domain artificially severs the pressure influence from outside the extraction window. Simply prescribing Dirichlet (extracted pressure values) or Neumann (zero gradient) conditions can lead to:
- Pressure drift over time
- Incorrect pressure gradients at the boundary
- Mass conservation issues
- Different flow behavior than the source simulation

The `spaceTimeWindowCoupledPressure` BC solves this by **blending Dirichlet and Neumann conditions** based on local flux. Two blending modes are available:

**fluxMagnitude mode (default, recommended):**
- **At stagnant faces (zero flux)**: Uses Dirichlet (safe to prescribe pressure)
- **At active flow faces**: Uses Neumann (gradient for continuity compatibility)

**flowDirection mode (original):**
- **At inflow faces**: More weight towards Dirichlet (pressure enters from outside)
- **At outflow faces**: More weight towards Neumann (pressure from interior dominates)

**Extraction configuration:**
```cpp
extractSubset
{
    type            spaceTimeWindowExtract;
    libs            (spaceTimeWindow);
    box             ((0.05 -0.15 0.01) (0.70 0.15 0.38));
    outputDir       "../subset";

    fields          (U p);          // Time-varying BC data (U and p for boundaries)
    initialFields   (U p nut);      // Initial conditions (nut needed for turbulence model)

    // Enable gradient extraction for pressure coupling
    extractGradients    true;
    gradientFields      (p);        // Extract normal gradient of pressure (creates gradp)

    writeFormat     deltaVarint;
    writeControl    timeStep;
    writeInterval   1;
}
```

**Reconstruction boundary conditions:**
```cpp
// In 0/p
oldInternalFaces
{
    type            spaceTimeWindowCoupledPressure;
    dataDir         "constant/boundaryData";
    phi             phi;
    blendingMode    fluxMagnitude;  // or flowDirection (default: fluxMagnitude)
    dirichletWeight 0.8;    // Weight at zero-flux faces (fluxMagnitude) or inflow (flowDirection)
    neumannWeight   0.8;    // Weight at active-flow faces (fluxMagnitude) or outflow (flowDirection)
    transitionWidth 0.1;    // Smooth transition width
    timeInterpolationScheme cubic;
    value           uniform 0;
}

// In 0/U
oldInternalFaces
{
    type            spaceTimeWindow;
    dataDir         "constant/boundaryData";
    timeInterpolationScheme cubic;
    value           uniform (0 0 0);
}
```

#### 3. Inlet-Outlet BC (`-inletOutletBC`) - For Quick Tests

Uses `spaceTimeWindowInletOutlet` BC which applies:
- **Velocity (U)**: Flux-based switching - Dirichlet (prescribed value) at inflow faces, zeroGradient at outflow faces
- **All scalar fields (p, nut, k, epsilon, omega)**: zeroGradient

This approach is simpler but **may produce different flow behavior** than the source simulation because pressure is not coupled to the exterior. Use for quick tests or when exact reproduction is not required.

```bash
spaceTimeWindowInitCase -sourceCase ../source-case -inletOutletBC
```

#### 4. Fixed Outlet Direction (`-outletDirection`) - For Steady-Mean Flows

Creates a separate `outlet` patch from faces at one edge of the extraction box:
- **oldInternalFaces**: Pure Dirichlet BC (spaceTimeWindow) on all remaining faces
- **outlet**: Pressure relief patch with inletOutlet for U and zeroGradient for p

Use when the mean flow direction is well-defined and outflow always occurs at a known location.

```bash
spaceTimeWindowInitCase -sourceCase ../source-case -outletDirection "(1 0 0)"
```

**Note:** Options `-preciseBC`, `-inletOutletBC`, and `-outletDirection` are mutually exclusive. Using incompatible combinations produces an error with guidance.

### Initial Fields vs Boundary Data Fields

The extraction can separate which fields are used for:
- **Initial conditions** (`initialFields`): Fields written to the start time directory for solver initialization
- **Time-varying boundary data** (`fields`): Fields written to `boundaryData/` for time-varying BCs

For precise BC (highest fidelity, requires `pimpleFoamPrecise`):

```cpp
// In extraction function object - precise BC (highest fidelity)
fields          (HbyA p);         // HbyA (intermediate velocity) and p
initialFields   (U p nut);        // U (final velocity) for solver initialization
extractGradients    true;
gradientFields      (p);          // Creates gradp for Neumann blending
```

For pressure-coupled BC (standard), extract both U and p with gradients:

```cpp
// In extraction function object - pressure-coupled (standard)
fields          (U p);             // Both U and p for time-varying BC
initialFields   (U p nut);         // nut needed for turbulence model
extractGradients    true;
gradientFields      (p);           // Creates gradp for Neumann blending
```

For quick tests with `-inletOutletBC`, only U needs time-varying data:

```cpp
// In extraction function object - quick tests only
fields          (U);               // Only U for time-varying BC (saves storage)
initialFields   (U p nut);         // More fields for initial conditions
```

Fields NOT in `boundaryData` automatically get `zeroGradient` BC on `oldInternalFaces`.

## Field Selection by Turbulence Model

When configuring the extraction, include all fields required by your turbulence model:

| Turbulence Model | Recommended Fields |
|------------------|-------------------|
| LES Smagorinsky | `(U p nut)` |
| LES dynamicKEqn | `(U p nut k)` |
| LES WALE | `(U p nut)` |
| RANS k-ε | `(U p nut k epsilon)` |
| RANS k-ω SST | `(U p nut k omega)` |
| Spalart-Allmaras | `(U p nut nuTilda)` |

For `-preciseBC` (highest fidelity), use `fields (HbyA p)` and `initialFields` as shown above. For standard pressure-coupled BC, use `fields (U p)` plus turbulence fields.

## Workflow

### Serial Execution

```bash
# 1. Run extraction during simulation
cd source-case
pimpleFoamPrecise    # Use pimpleFoamPrecise for -preciseBC, or pimpleFoam otherwise

# 2. Initialize reconstruction case
cd ../subset-case
spaceTimeWindowInitCase -sourceCase ../source-case -preciseBC  # or no flag for standard

# 3. Run reconstruction
pimpleFoam
```

In serial mode:
- Mesh and initial fields are written directly to the output directory
- Boundary data is written at every timestep
- The extractor forces field writes at t_1 and t_2 (for interpolation buffers)

### Parallel Execution

Both the source simulation (extraction phase) and the reconstruction simulation can run in parallel. This enables efficient use of HPC resources for both phases of the workflow.

```bash
# 1. Run extraction during parallel simulation
cd source-case
mpirun -np 8 pimpleFoam -parallel    # Extraction happens automatically

# 2. Reconstruct the start time (t_2 for cubic interpolation)
reconstructPar -time 0.0002

# 3. Initialize reconstruction case (recommended: inlet-outlet BC)
cd ../subset-case
spaceTimeWindowInitCase -sourceCase ../source-case -inletOutletBC

# 4. Run reconstruction (serial or parallel)
pimpleFoam                           # Serial
# OR
decomposePar && mpirun -np 4 pimpleFoam -parallel   # Parallel
```

In parallel mode:
- Boundary data from all processors is gathered to master and written as single files
- The extractor forces field writes at t_0, t_1, and t_2 for interpolation buffers (t_0 only in parallel)
- Only extraction box parameters are written (not mesh) - `spaceTimeWindowInitCase` creates the mesh
- `reconstructPar` must be run at the reconstruction start time (t_2) to provide source fields

**Tip:** The reconstruction case is typically much smaller than the source case, so fewer processors may be needed. The spaceTimeWindow boundary conditions work identically in serial and parallel modes.

### Time Interpolation Modes

| Mode | Start Time | Timesteps Used | Use Case |
|------|------------|----------------|----------|
| `none` (exact) | t_0 | 1 (exact match) | Highest fidelity (no interpolation error) |
| `linear` | t_1 | 2 (bracketing) | Simple, smoothly varying flows |
| `cubic` | t_2 | 4 (Catmull-Rom) | Unsteady turbulent flows (recommended) |

`spaceTimeWindowInitCase` automatically sets `startTime` and `endTime` to ensure sufficient buffer timesteps for the selected interpolation scheme.

## Components

### spaceTimeWindowExtract (Function Object)

Extracts face-interpolated boundary values at every timestep during the original simulation.

```cpp
functions
{
    extractSubset
    {
        type            spaceTimeWindowExtract;
        libs            (spaceTimeWindow);

        // Define extraction region with bounding box (min max points)
        // Box must be fully internal to domain (no intersection with boundaries)
        box             ((0.05 -0.25 0.01) (0.90 0.25 0.38));

        outputDir       "../subset-case";        // Output case directory

        // Fields for time-varying boundary data (written every timestep)
        fields          (U p);                   // U and p for pressure-coupled BC

        // Fields for initial conditions (optional, defaults to 'fields')
        initialFields   (U p nut);               // nut needed for turbulence model

        // Enable gradient extraction for pressure-coupled BC
        extractGradients    true;
        gradientFields      (p);                 // Creates gradp for Neumann part

        // Write format for boundary data (optional)
        // Options: ascii, binary, deltaVarint, dvzt
        // Default: ascii
        writeFormat     dvzt;

        // Precision for deltaVarint (optional, default: 6)
        deltaVarintPrecision  6;

        // Compression for ascii/binary (optional, ignored for deltaVarint)
        // writeCompression on;

        writeControl    timeStep;
        writeInterval   1;
    }
}
```

**Parameters:**

| Parameter    | Description                                    | Required | Default |
|--------------|------------------------------------------------|----------|---------|
| box          | Bounding box as ((minX minY minZ) (maxX maxY maxZ)) | yes | - |
| outputDir    | Output case directory for extracted data       | yes      | - |
| fields       | List of fields for time-varying BC (boundaryData) | yes   | - |
| initialFields | List of fields for initial conditions         | no       | same as fields |
| writeFormat  | Format for boundary data: `ascii`, `binary`, `deltaVarint`, or `dvzt` | no | ascii |
| deltaVarintPrecision | Decimal digits for delta-varint quantization | no | 6 |
| writeCompression | Gzip compression for ascii/binary (ignored for deltaVarint) | no | off |

**Field Selection Strategy:**

For accurate reconstruction with pressure-coupled BC (recommended):

```cpp
fields          (U p);             // Both U and p for time-varying BC
initialFields   (U p nut);         // All fields needed for solver startup
extractGradients    true;
gradientFields      (p);           // Extract normal gradient for pressure coupling
```

For quick tests with `-inletOutletBC` (not recommended for accurate reconstruction):

```cpp
fields          (U);               // Only velocity needs time-varying BC data
initialFields   (U p nut);         // All fields needed for solver startup
```

**Output structure:**
```
outputDir/
    constant/
        polyMesh/               # Subset mesh (only oldInternalFaces patch)
        boundaryData/
            oldInternalFaces/
                points              # Face centres (reference only)
                extractionMetadata  # Settings and timestep list
                <time>/
                    U               # or U.dvz, U.dvz.zstd, U.dvzt.zstd, U.dvz.enc, U.dvz.zstd.enc
    <startTime>/                # Initial subset fields (always ASCII)
        U
        p
        nut
```

### spaceTimeWindowInitCase (Utility)

Initializes a fully configured reconstruction case from the extracted data.

```bash
# Highest fidelity: precise BC (requires pimpleFoamPrecise extraction)
spaceTimeWindowInitCase -sourceCase ../source-case -preciseBC

# Standard: pressure-coupled BC (default, no flag needed)
spaceTimeWindowInitCase -sourceCase ../source-case

# Quick tests: inlet-outlet BC
spaceTimeWindowInitCase -sourceCase ../source-case -inletOutletBC

# Alternative: fixed outlet direction for steady-mean flows
spaceTimeWindowInitCase -sourceCase ../source-case -outletDirection "(1 0 0)"

# With mass flux correction (optional, ensures exact mass conservation)
spaceTimeWindowInitCase -sourceCase ../source-case -inletOutletBC -correctMassFlux
```

**Options:**

| Option           | Description                                      | Required |
|------------------|--------------------------------------------------|----------|
| -sourceCase      | Source case directory (where extraction ran)     | yes      |
| -extractDir      | Directory with extracted data (default: cwd)     | no       |
| -preciseBC       | Use HbyA-based velocity BC and coupledPressure for p (requires `pimpleFoamPrecise` extraction) | no |
| -inletOutletBC   | Use flux-based inlet-outlet BC for U, zeroGradient for scalars (for quick tests) | no |
| -outletDirection | Create fixed outlet patch in given direction (e.g., "(1 0 0)") | no |
| -outletFraction  | Fraction of box extent for outlet region (default: 0.1) | no |
| -correctMassFlux | Apply least-squares mass flux correction to boundaryData | no |
| -initialFields   | Override initial fields list (e.g., "(U p nut k)") | no |
| -refineLevel     | Refine mesh N times (each level splits cells ×8) | no |
| -coarsenLevel    | Coarsen mesh N times (each level merges cells)   | no |
| -overwrite       | Overwrite existing files                         | no       |

**Note:** `-preciseBC`, `-inletOutletBC`, and `-outletDirection` are mutually exclusive. `-refineLevel` and `-coarsenLevel` are also mutually exclusive.

**What it creates:**

1. `system/controlDict` - With matching solver, deltaT, adjustTimeStep from extraction
   - `startTime` set to t_2 (third timestep) for cubic interpolation buffer
   - `endTime` set to t_{n-2} (third-to-last) for cubic interpolation buffer
2. `system/fvSchemes`, `system/fvSolution` - Copied from source case (with pRefPoint added)
3. `constant/` files - All physics properties copied (mandatory for fidelity)
4. Initial field files with appropriate BCs based on options and boundaryData availability

**Boundary Condition Assignment Logic:**

For each field in the initial time directory:

| Mode | Field | In boundaryData? | Field Type | BC Applied |
|------|-------|-------------------|------------|------------|
| `-preciseBC` | U | --- | vector | `spaceTimeWindow` + `fieldTableName HbyA` |
| `-preciseBC` | p | --- | scalar | `spaceTimeWindowCoupledPressure` |
| `-preciseBC` | other | No | any | `zeroGradient` |
| `-inletOutletBC` | --- | Yes | vector (U) | `spaceTimeWindowInletOutlet` |
| `-inletOutletBC` | --- | Yes | scalar | `zeroGradient` |
| default | --- | Yes | any | `spaceTimeWindow` (Dirichlet) |
| any | --- | No | any | `zeroGradient` |

This ensures fields without time-varying data automatically get appropriate BCs. With `-preciseBC`, the U and p fields are handled specially (U reads HbyA data via `fieldTableName`).

### pimpleFoamPrecise (Solver)

Modified version of OpenFOAM's `pimpleFoam` that registers the intermediate velocity **u\*** = H(U)/A (named `HbyA`) in the objectRegistry, making it accessible to function objects like `spaceTimeWindowExtract`.

**Why a separate solver is needed:**

In standard `pimpleFoam`, the local `volVectorField HbyA(...)` in `pEqn.H` is stack-allocated. It auto-registers in the objectRegistry upon construction but its destructor **unregisters** it when `pEqn.H` goes out of scope. By the time function objects run (after `runTime.write()`), HbyA is no longer available.

**Implementation:**

`pimpleFoamPrecise` solves this by declaring a persistent `volVectorField HbyA_stored` in `createFields.H` with `IOobject::REGISTER`. This field survives across the time loop. In `pEqn.H`, the computed HbyA is copied to the persistent field with `HbyA_stored = HbyA`.

The overhead is negligible: one field copy per pressure correction iteration.

**Usage:**

```bash
# In source case controlDict:
#   application  pimpleFoamPrecise;
# In extraction function object:
#   fields  (HbyA p);

cd source-case
pimpleFoamPrecise
```

The reconstruction phase uses standard `pimpleFoam` --- no special solver is needed for reconstruction.

### spaceTimeWindowCoupledPressure (Boundary Condition)

Mixed boundary condition that blends Dirichlet (extracted pressure) and Neumann (extracted gradient) based on local flux direction. This is the **recommended approach** for pressure because it properly handles the elliptic nature of the pressure equation.

```cpp
oldInternalFaces
{
    type            spaceTimeWindowCoupledPressure;
    dataDir         "constant/boundaryData";
    phi             phi;                        // Flux field name
    blendingMode    fluxMagnitude;              // or flowDirection
    dirichletWeight 0.8;                        // Dirichlet strength (0-1)
    neumannWeight   0.8;                        // Neumann strength (0-1)
    transitionWidth 0.1;                        // Smooth transition width
    allowTimeInterpolation  true;
    timeInterpolationScheme cubic;              // or "linear" or "none"
    value           uniform 0;
}
```

**How it works:**

The BC implements a mixed (Robin) condition that adapts based on local flux:

1. At each timestep, reads both pressure values (`p`) and normal gradients (`gradp`) from boundaryData
2. Computes blending weight `w` for each face based on flux using a smooth tanh transition
3. Applies: `p_face = w * p_Dirichlet + (1-w) * (p_cell + gradp_Neumann * d)`

**fluxMagnitude mode** (recommended):
- At zero-flux faces: `w -> dirichletWeight` (safe to prescribe pressure at stagnant faces)
- At high-flux faces: `w -> 1 - neumannWeight` (use gradient for continuity compatibility)

**flowDirection mode** (original):
- At inflow faces (phi < 0): `w -> dirichletWeight` (pressure enters from outside)
- At outflow faces (phi > 0): `w -> 1 - neumannWeight` (pressure from interior dominates)

**Mathematical formulation:**

The blending weight uses a smooth hyperbolic tangent transition:
```
transition = tanh(phi_normalized / transitionWidth)
w = 0.5 * ((w_dir + w_neu) - (w_dir - w_neu) * transition)
```

This ensures C∞ continuity at the transition, avoiding numerical artifacts.

**Properties:**

| Property       | Description                                    | Required | Default              |
|----------------|------------------------------------------------|----------|----------------------|
| dataDir        | Path to boundaryData directory                 | no       | constant/boundaryData |
| phi            | Name of flux field                             | no       | phi                  |
| blendingMode   | `fluxMagnitude` or `flowDirection`             | no       | fluxMagnitude        |
| dirichletWeight | Dirichlet contribution strength (0-1)         | no       | 0.8                  |
| neumannWeight  | Neumann contribution strength (0-1)            | no       | 0.8                  |
| transitionWidth | Normalized flux width for smooth transition   | no       | 0.1                  |
| pressureField  | Name of pressure data file in boundaryData     | no       | p                    |
| gradientField  | Name of gradient data file in boundaryData     | no       | gradp                |
| allowTimeInterpolation | Permit interpolation for missing timesteps | no | true |
| timeInterpolationScheme | `none`, `linear`, or `cubic`            | no       | cubic                |

**Extraction requirements:**

To use this BC, extraction must include gradient data:
```cpp
extractSubset
{
    // ... other settings ...
    fields              (U p);
    extractGradients    true;
    gradientFields      (p);    // Extracts "gradp" (normal gradient)
}
```

### spaceTimeWindowInletOutlet (Boundary Condition)

Flux-based boundary condition that reads pre-computed velocity values and applies them only at inflow faces. Use for velocity when using the pressure-coupled approach, or for quick tests.

```cpp
oldInternalFaces
{
    type            spaceTimeWindowInletOutlet;
    dataDir         "constant/boundaryData";
    phi             phi;                        // Flux field name
    allowTimeInterpolation  true;
    timeInterpolationScheme cubic;              // or "linear" or "none"
    value           uniform (0 0 0);
}
```

**How it works:**
1. At each timestep, reads velocity values from boundaryData (with optional time interpolation)
2. Computes flux through each face: `phi_f = U_f . Sf`
3. For inflow faces (phi < 0): applies prescribed velocity from boundaryData
4. For outflow faces (phi >= 0): applies zeroGradient (extrapolates from interior)

**Properties:**

| Property       | Description                                    | Required | Default              |
|----------------|------------------------------------------------|----------|----------------------|
| dataDir        | Path to boundaryData directory                 | no       | constant/boundaryData |
| phi            | Name of flux field                             | no       | phi                  |
| allowTimeInterpolation | Permit interpolation for missing timesteps | no | false |
| timeInterpolationScheme | `none`, `linear`, or `cubic`            | no       | linear               |

### spaceTimeWindow (Boundary Condition)

Pure Dirichlet boundary condition that prescribes values on all faces regardless of flow direction.

```cpp
oldInternalFaces
{
    type            spaceTimeWindow;
    dataDir         "constant/boundaryData";
    fixesValue      true;                       // Tells adjustPhi these values are fixed
    allowTimeInterpolation  true;
    timeInterpolationScheme cubic;
    value           uniform (0 0 0);
}
```

Use this for:
- Scalar fields when not using `-inletOutletBC`
- Situations where fixed values are explicitly desired on all faces

**Properties:**

| Property       | Description                                    | Required | Default              |
|----------------|------------------------------------------------|----------|----------------------|
| dataDir        | Path to boundaryData directory                 | no       | constant/boundaryData |
| fixesValue     | Report to adjustPhi that values are fixed      | no       | true                 |
| allowTimeInterpolation | Permit interpolation for missing timesteps | no | false |
| timeInterpolationScheme | `none`, `linear`, or `cubic`            | no       | linear               |
| reportFlux     | Print net flux through patch (velocity only)   | no       | false                |

### Time Interpolation Schemes

| Scheme | Description | Points Used | Buffer Needed |
|--------|-------------|-------------|---------------|
| `none` | Exact timestep matching required | 1 | 0 |
| `linear` | Linear interpolation between bracketing timesteps | 2 | 1 at each end |
| `cubic` | Centripetal Catmull-Rom cubic spline | 4 | 2 at each end |

**Cubic interpolation** is recommended for turbulent flows because:
- Provides C1 continuity (smooth first derivatives)
- Uses centripetal parameterization to handle non-uniform time spacing correctly
- Essential for adaptive timestepping (`adjustTimeStep = yes`)
- Better captures temporal variations without overshoots

## Mass Conservation

### The Challenge

Face interpolation during extraction can introduce small mass imbalances. Additionally, with all-Dirichlet BCs, `adjustPhi()` cannot correct these imbalances.

### Solutions

**1. Use Pressure-Coupled BC (Recommended)**

The `spaceTimeWindowCoupledPressure` BC combined with `spaceTimeWindow` for velocity naturally handles mass conservation through proper pressure coupling.

**2. Use `-inletOutletBC` (For Quick Tests)**

The flux-based inlet-outlet BC handles mass conservation:
- Outflow faces use zeroGradient, allowing natural outflow
- No artificial mass imbalance from prescribed outflow velocities
- Works without any special mass correction

**3. Use `-correctMassFlux`**

Applies least-squares correction to boundaryData to ensure exact mass conservation:

```bash
spaceTimeWindowInitCase -sourceCase ../source -inletOutletBC -correctMassFlux
```

The correction minimizes ||U_corrected - U||² subject to Σ(U_corrected · Sf) = 0:
```
U_corrected = U - (imbalance / totalSfMag) * n
```

**4. Use `-outletDirection` with `fixesValue true`**

Creates an outlet patch where mass imbalance can escape:
- `oldInternalFaces` uses `fixesValue true` (values not modified by adjustPhi)
- `outlet` uses `inletOutlet` BC (allows adjustPhi correction)

### Recommended Configurations

| Use Case | Command/Approach | Notes |
|----------|---------|-------|
| Highest fidelity reconstruction | `-preciseBC` | ~20% lower pressure error than default. Requires `pimpleFoamPrecise`. |
| Accurate reconstruction | Pressure-coupled BC (default) | Reproduces original flow to solver tolerance. |
| Quick tests | `-inletOutletBC` | Natural mass balance, may drift from source. |
| Steady-mean flow direction | `-outletDirection "(1 0 0)"` | When outlet location is known. |
| Maximum fidelity | `-preciseBC` + `-correctMassFlux` | Best accuracy. |

## Boundary Data Compression

The `writeFormat` parameter controls how boundary data files are written.

| Format | Extension | Typical Size | Notes |
|--------|-----------|--------------|-------|
| `ascii` | (none) | 100% baseline | Human-readable, default |
| `binary` | (none) | ~50% | OpenFOAM native binary |
| `ascii` + gzip | .gz | ~10% | `writeCompression on` |
| `binary` + gzip | .gz | ~8% | Combined |
| `deltaVarint` | .dvz | **~2.7%** | High compression, self-contained |
| `deltaVarint` + zstd | .dvz.zstd | **~0.4%** | Integrated zstd (auto if libzstd present) |
| `dvzt` | .dvzt | **~2.4%** | Best compression, recommended |
| `dvzt` + zstd | .dvzt.zstd | **~0.3%** | Best compression + integrated zstd |

### Delta-Varint Codec (DVZ)

Specialized codec optimized for CFD boundary data:

1. Component-major ordering (groups similar values)
2. **Spatial delta encoding** (stores differences between consecutive face values within the same timestep)
3. Quantization to configurable precision
4. Variable-length integer encoding (zigzag encoding for signed deltas)

**Important:** Each timestep file is completely self-contained. Delta encoding is purely spatial (between consecutive faces in a single field), not temporal. This means:
- Any timestep can be read independently without access to other timesteps
- No temporal decompression dependencies
- Safe to delete individual timestep files without affecting others

```cpp
writeFormat     deltaVarint;
deltaVarintPrecision  6;    // ~1e-6 relative precision
```

### Integrated Zstd Compression (Optional)

When built with libzstd support (auto-detected by `./Allwmake`), DVZ and DVZT data is automatically wrapped with zstd level-3 compression before writing to disk. This is applied transparently as a compression layer:

```
Encoding pipeline:  codec encode → zstd compress → [optional encrypt] → write
Extension stacking: fieldName.dvz.zstd  or  fieldName.dvzt.zstd  or  fieldName.dvz.zstd.enc
```

The zstd layer provides an additional ~10x reduction on top of the already-compact varint encoding, bringing typical file sizes to ~0.3-0.4% of uncompressed ASCII. Decompression is transparent at runtime - the boundary conditions automatically detect and decompress `.zstd` files.

**No configuration needed**: zstd compression is applied automatically when the library is built with libzstd. Files without `.zstd` extension (from older builds) remain fully readable for backward compatibility.

### Delta-Varint-Temporal Codec (DVZT)

Enhanced codec that exploits both spatial and temporal correlation for better compression:

1. **Keyframes** (every N timesteps): Self-contained, same as DVZ
2. **Delta frames**: Hybrid spatial-temporal prediction
   - Uses weighted prediction: `predicted = 0.3 * spatial_neighbor + 0.7 * temporal_neighbor`
   - Encodes residuals (actual - predicted) instead of raw spatial deltas
   - Typically ~10% smaller than DVZ for delta frames

```cpp
writeFormat     dvzt;
deltaVarintPrecision  6;        // ~1e-6 relative precision
dvztKeyframeInterval  20;       // Keyframe every 20 timesteps (default)
```

**DVZT Workflow:**

During extraction, DVZT writes smaller .dvzt files (or .dvzt.zstd with integrated zstd). During case initialization, `spaceTimeWindowInitCase` automatically converts .dvzt files to .dvz format (required because delta frames need sequential processing). The zstd layer is preserved through the conversion. The resulting .dvz (or .dvz.zstd) files are read by the spaceTimeWindow BC at runtime.

```
Extraction:  .dvzt.zstd files (with zstd) or .dvzt files (without)
     |
     v
spaceTimeWindowInitCase:  Converts .dvzt.zstd -> .dvz.zstd  (or .dvzt -> .dvz)
     |
     v
Runtime:  Reads .dvz.zstd or .dvz files (auto-detected by BC)
```

**When to use DVZT:**
- Long simulations with many timesteps (storage savings accumulate)
- Network/disk bandwidth constraints during extraction
- Small timesteps (Δt = 1e-5 or 1e-6) where temporal correlation is very strong
- DNS or acoustic simulations requiring fine temporal resolution

**Compression vs Timestep Size:**

DVZT benefits increase dramatically with smaller timesteps because consecutive values become nearly identical:

| Δt | Temporal Correlation | DVZT vs DVZ Savings |
|----|---------------------|---------------------|
| 1e-4 | Moderate | ~10% smaller |
| 1e-5 | Strong | ~20-30% smaller |
| 1e-6 | Very strong | ~30-50% smaller |

With very small Δt, most temporal residuals quantize to near-zero values, requiring only 1 byte each in varint encoding.

**When to use DVZ:**
- Simpler workflow (no conversion step)
- When random access to individual timesteps is needed during extraction
- Shorter simulations where DVZT overhead isn't worth it
- Large timesteps where temporal correlation is weak

## Encryption (Optional)

Boundary data can be encrypted using X25519 asymmetric encryption (libsodium sealed boxes). When both zstd and encryption are enabled, the extension stacking is `fieldName.dvz.zstd.enc` (or `.dvzt.zstd.enc`). During case initialization, decryption strips the `.enc` extension, leaving the zstd-compressed file for the BC to read.

### Building with Encryption Support

```bash
export FOAM_USE_SODIUM=1 && ./Allwmake
```

### Usage

```bash
# Generate key pair
spaceTimeWindowKeygen
# Public key:  fqzYQ0U8j27tFEr5WzEMylbvXYP+9CAyk0JhwwZ2rwg=
# Private key: ************H8exb*****n4Kv0nu5gqljI2RPBCwl4=

# Extraction with public key
# Add to function object: publicKey "fqzYQ0U8j27tFEr5WzEMylbvXYP+9CAyk0JhwwZ2rwg=";

# Decryption during case init (prompts for private key)
spaceTimeWindowInitCase -sourceCase ../source-case -inletOutletBC
# Enter private key (base64):
```

### Secure Key Storage

The private key should be stored securely:
- **Hardware security keys**: YubiKey, Nitrokey, or similar FIDO2/PIV-capable devices (strongest protection)
- **Secure password managers**: KeePassXC, 1Password, Bitwarden with strong master passwords
- **Encrypted vaults**: GPG-encrypted files or OS keychains (macOS Keychain, GNOME Keyring)

**Warning**: Never commit the private key to version control.

### Ethical Considerations and HPC Transparency

Most high-performance computing (HPC) centers require full transparency regarding simulations performed on their infrastructure. For this reason, **encryption support is an optional compile-time feature** controlled by the `FOAM_USE_SODIUM` flag. Centers may choose not to enable this functionality.

**What encryption protects:**
- The time-varying boundary field data stored in `constant/boundaryData`. When downloaded to a client workstation, this information is sealed and cannot be read without the private key.

**What encryption does NOT protect:**
- The test case setup (mesh, dictionaries, initial conditions) remains unencrypted. Anyone can re-run the global simulation to regenerate the boundary data - encryption does not prevent this.

**Practical security model:**

1. Re-running the global simulation to obtain internal fields is **computationally expensive** and requires HPC resources. Such jobs are recorded in scheduler logs, providing an audit trail.
2. Offline recalculation on a personal workstation is typically infeasible due to computational cost.
3. If file permissions on HPC scratch directories are misconfigured and the user chooses not to write timestep data from the global simulation, other users can only copy the unencrypted test case setup - not the valuable boundary data.

**Note**: The encryption feature is designed for protecting intellectual property during data distribution, not for hiding simulations from HPC administrators. Always comply with your computing center's acceptable use policies.

### Protection Against Malware and Ransomware

The spaceTimeWindow approach offers an additional security benefit: by storing only the minimal data required for reconstruction, the valuable simulation results can be protected against malware and ransomware attacks.

**Secure archival strategy:**

The reconstruction case contains only:
- The subset mesh (`constant/polyMesh`)
- Encrypted boundary data (`constant/boundaryData/*.enc`)
- Initial fields for one timestep
- Configuration dictionaries

This minimal dataset (typically a few gigabytes) can be stored on:
- **Write-Once Read-Many (WORM) media**: M-DISC, tape archives, or cloud storage with object lock (AWS S3 Object Lock, Azure Immutable Blob Storage)
- **Secure vaults**: Hardware-encrypted drives kept offline
- **Air-gapped backups**: Disconnected storage unreachable by network-based malware

**Security model:**
1. Encrypted boundary data on WORM storage cannot be modified or deleted by malware
2. Private key stored separately (hardware token or password manager), not on the same system
3. Reconstruction always possible - archived data remains intact even if systems are compromised
4. Small archive size makes backup and verification practical, unlike terabyte-scale full simulation outputs

**Tip**: For critical simulations, archive the encrypted reconstruction case to WORM storage immediately after extraction. The full simulation data can then be deleted from HPC scratch space.

## Building

```bash
# Build everything (library + utilities)
# Auto-detects libzstd (compression) and libsodium (encryption)
./Allwmake
```

Optional dependencies are auto-detected by `./Allwmake`:
- **libzstd**: Enables integrated zstd compression of DVZ/DVZT files (~10x additional reduction)
- **libsodium**: Enables X25519 encryption of boundary data

### Manual Build

```bash
# Build the library
cd src/spaceTimeWindow
wmake

# Build utilities
cd applications/utilities/preProcessing/spaceTimeWindowInitCase
wmake
```

For manual builds with optional features, export the flags before `wmake`:
```bash
export FOAM_USE_ZSTD=1    # if libzstd is installed
export FOAM_USE_SODIUM=1  # if libsodium is installed
```

## Documentation

API documentation can be generated with Doxygen:

```bash
./Allwmake doc
```

## Requirements

- **OpenFOAM v2512** (openfoam.com) or compatible ESI-OpenCFD version
- For parallel extraction: MPI environment
- For integrated zstd compression: libzstd development package (libzstd-dev or libzstd-devel)
- For encryption: libsodium development package (libsodium-dev or libsodium-devel)

## Example Case

The `examples/ufr2-02` directory contains a complete LES test case based on the ERCOFTAC Classic Collection Database Case 043:

- **Case**: ERCOFTAC UFR2-02 - Flow around a square cylinder at Re = 21,400
- **Reference**: Lyn et al. (1995) J. Fluid Mech. 304, 285-319
- **Mesh**: Based on the mesh generation script by Niklas Nordin from the [OpenFOAM Wiki UFR2-02 Benchmark](https://openfoamwiki.net/index.php?title=Benchmark_ercoftac_ufr2-02)
- **ERCOFTAC Database**: http://cfd.mace.manchester.ac.uk/ercoftac/ (Case 043)

The case demonstrates vortex shedding from a square cylinder with a von Karman vortex street in the wake region - an ideal scenario for space-time window extraction since the extraction box can capture the wake dynamics while excluding the cylinder and inlet regions.

```bash
cd examples/ufr2-02
./Allrun
```

## Validation

### CFD Validation Philosophy

In computational fluid dynamics, validation means demonstrating that a simulation reproduces the **physics** of the flow—not that it produces bit-identical numbers. A validated CFD result agrees with experimental measurements (velocity profiles, pressure distributions, turbulent statistics, shedding frequencies, etc.) within quantified uncertainties.

This distinction is central to the space-time window approach. The reconstruction does not aim to be a bit-exact copy of the original simulation; it aims to be a **physically faithful re-simulation** of the same flow inside the region of interest. The governing equations, boundary conditions, and fluid properties are identical; what changes is the domain extent and, optionally, the mesh resolution.

Different meshes, turbulence models, numerical schemes, and even different software packages will produce numerically different fields, yet all can be equally "valid" if they capture the physical phenomena to the required accuracy.

### Validating the Reconstruction

The reconstruction can be validated by comparing:

- **Instantaneous velocity fields** at matching timesteps
- **Vortex shedding frequency** (should match original)
- **Mean velocity profiles** in the wake
- **Reynolds stress components**
- **Pressure distributions** and recovery

Because the reconstruction solves the same equations with the same physical parameters, any discrepancy with the full-domain solution is due to the finite-volume discretization change at the extraction boundary (see Accuracy Comparison above), not to different physics. These errors are localized near the boundary and decay rapidly with distance.

**Note:** The reconstruction should closely reproduce the original flow within the extraction region (to solver tolerance) when using exact timestep matching. Small differences may occur at the outlet boundary due to the different boundary condition type. When the mesh resolution is changed (refinement or coarsening), additional spatial interpolation errors appear, but the physical behavior—vortex dynamics, mean profiles, spectral content—is preserved.

## Mesh Coarsening and Refinement

The reconstruction mesh can be coarsened or refined relative to the extraction mesh, enabling simulations at different resolutions than the original.

```bash
# Refine mesh (finer than extraction)
spaceTimeWindowInitCase -sourceCase ../source -inletOutletBC -refineLevel 1

# Coarsen mesh (coarser than extraction)
spaceTimeWindowInitCase -sourceCase ../source -inletOutletBC -coarsenLevel 1

# Multiple levels (each level = 8× cell count change)
spaceTimeWindowInitCase -sourceCase ../source -inletOutletBC -refineLevel 2
```

### Spatial Interpolation Algorithms

When the reconstruction mesh differs from the extraction mesh, spatial interpolation is required for boundary data. The spaceTimeWindow boundary conditions handle this automatically at runtime:

**Refinement** (more target faces than source): **Barycentric interpolation with 2D Delaunay triangulation**
- Source face centers are triangulated using the Bowyer-Watson algorithm
- For each target face center, the enclosing triangle is found and barycentric weights computed
- If the target point lies outside all triangles (due to irregular submesh boundaries from original cell shapes), the algorithm finds the nearest triangle by centroid distance and uses clamped barycentric coordinates
- Provides smooth C⁰ continuous interpolation

**Coarsening** (fewer target faces than source): **Area-weighted averaging**
- An octree is built from source points for efficient spatial lookup
- For each target face, all source points within a search radius are averaged with equal weights
- Search radius is progressively expanded if no points are found
- Ensures conservation of integral quantities

The interpolation works on 2D point clouds. Points are grouped by which face of the extraction bounding box they belong to (6 planar surfaces: ±X, ±Y, ±Z), then projected to 2D by dropping the constant coordinate. This handles the irregular face distribution that arises from extracting a submesh with original cell shapes.

### Initial Field Interpolation

Initial fields are also interpolated when mesh resolution changes:
- **Refinement**: Uses `mapFields` with cell-center interpolation
- **Coarsening**: Uses volume-weighted averaging of source cells

**Note:** Spatial interpolation introduces smoothing, particularly for coarsening. For turbulent flows, this may affect small-scale structures. Consider the trade-off between computational cost and resolution fidelity.

## External Compression Benchmark

> **Note:** With integrated zstd compression (built with libzstd), boundary data is automatically compressed with zstd -3 on the fly. The benchmarks below show what is achievable with *external* tools on *uncompressed* DVZ/DVZT files. When integrated zstd is enabled, the files are already zstd-compressed, so applying external zstd again provides no benefit. However, external tools like 7z LZMA2 can still achieve better ratios for archival.

DVZ and DVZT files can be further compressed using external tools for archival or transfer. The following benchmarks compare various compression algorithms on real boundary data files.

### Test Environment

- **CPU**: Intel Core i7-4790 @ 3.60 GHz (4 cores, 8 threads, Haswell)
- **L3 Cache**: 8 MB
- **RAM**: DDR3-1600

### Test Data

| Format | Files | Total Size | Description |
|--------|-------|------------|-------------|
| DVZ | 618 | 13 MB | Spatial delta-varint encoded |
| DVZT | 938 | 30 MB | Temporal delta-varint encoded |

The compression results for DVZ files are shown in the DVZ Files Compression Results table below, and DVZT results in the DVZT Files Compression Results table.

### DVZ Files Compression Results

| Method | Size | Ratio | Speed | Notes |
|--------|------|-------|-------|-------|
| 7z lzma2 -mx9 | 769 KB | 6.01% | 11.2 MB/s | **Best ratio** |
| 7z lzma -mx9 | 770 KB | 6.02% | 11.4 MB/s | |
| xz -9 | 769 KB | 6.01% | 10.1 MB/s | |
| xz -6 | 769 KB | 6.01% | 10.2 MB/s | |
| xz -1 | 1.5 MB | 11.37% | 36.0 MB/s | |
| zstd --ultra -22 | 964 KB | 7.53% | 4.7 MB/s | Very slow |
| zstd -19 | 965 KB | 7.54% | 15.3 MB/s | |
| rar -m5 | 1017 KB | 7.95% | 26.9 MB/s | |
| rar -m3 | 1019 KB | 7.96% | 45.1 MB/s | |
| 7z ppmd -mx9 | 1.4 MB | 11.02% | 10.7 MB/s | |
| 7z lzma2 -mx1 | 1.5 MB | 11.79% | 137.6 MB/s | |
| 7z lzma -mx1 | 1.5 MB | 11.60% | 48.8 MB/s | |
| zstd -9 | 1.7 MB | 13.17% | 109.5 MB/s | |
| 7z ppmd -mx5 | 1.7 MB | 13.13% | 15.2 MB/s | |
| zstd -3 | 1.8 MB | 13.65% | 309.0 MB/s | **Best balance** |
| zstd -1 | 1.9 MB | 14.43% | 533.7 MB/s | |
| bzip2 -9 | 2.0 MB | 15.28% | 13.5 MB/s | |
| bzip2 -5 | 2.0 MB | 15.50% | 13.7 MB/s | |
| lz4 -9 | 2.1 MB | 16.11% | 152.7 MB/s | |
| gzip -9 | 2.1 MB | 16.38% | 62.9 MB/s | |
| gzip -6 | 2.1 MB | 16.39% | 94.2 MB/s | |
| gzip -1 | 2.2 MB | 17.17% | 146.4 MB/s | |
| bzip2 -1 | 2.2 MB | 17.57% | 13.9 MB/s | |
| lz4 | 2.3 MB | 17.80% | 635.2 MB/s | **Fastest** |

### DVZT Files Compression Results

| Method | Size | Ratio | Speed | Notes |
|--------|------|-------|-------|-------|
| 7z lzma2 -mx9 | 1.9 MB | 6.31% | 10.1 MB/s | **Best ratio** |
| 7z lzma -mx9 | 1.9 MB | 6.31% | 11.5 MB/s | |
| xz -9 | 1.9 MB | 6.31% | 9.3 MB/s | |
| xz -6 | 1.9 MB | 6.32% | 10.5 MB/s | |
| 7z lzma2 -mx5 | 2.0 MB | 6.40% | 17.1 MB/s | |
| 7z lzma -mx5 | 2.0 MB | 6.40% | 17.2 MB/s | |
| xz -1 | 2.4 MB | 7.97% | 42.3 MB/s | |
| 7z lzma -mx1 | 2.4 MB | 7.98% | 65.4 MB/s | |
| 7z lzma2 -mx1 | 2.5 MB | 8.24% | 219.4 MB/s | |
| zstd --ultra -22 | 2.6 MB | 8.48% | 2.4 MB/s | Very slow |
| zstd -19 | 2.6 MB | 8.50% | 10.4 MB/s | |
| rar -m5 | 2.8 MB | 9.28% | 28.0 MB/s | |
| rar -m3 | 2.8 MB | 9.29% | 46.8 MB/s | |
| zstd -9 | 2.8 MB | 9.27% | 120.6 MB/s | |
| 7z ppmd -mx9 | 2.8 MB | 9.37% | 12.2 MB/s | |
| zstd -3 | 3.1 MB | 10.10% | 396.1 MB/s | **Best balance** |
| zstd -1 | 3.1 MB | 10.25% | 614.0 MB/s | |
| 7z ppmd -mx5 | 3.1 MB | 10.10% | 20.6 MB/s | |
| bzip2 -9 | 3.5 MB | 11.52% | 13.7 MB/s | |
| lz4 -9 | 3.5 MB | 11.50% | 190.0 MB/s | |
| bzip2 -5 | 3.5 MB | 11.76% | 14.2 MB/s | |
| gzip -9 | 3.8 MB | 12.69% | 79.5 MB/s | |
| gzip -6 | 3.8 MB | 12.71% | 104.3 MB/s | |
| bzip2 -1 | 4.1 MB | 13.60% | 13.9 MB/s | |
| lz4 | 4.1 MB | 13.56% | 763.0 MB/s | **Fastest** |
| gzip -1 | 4.4 MB | 14.50% | 163.7 MB/s | |

### Recommendations

| Use Case | Method | Ratio | Speed | Command |
|----------|--------|-------|-------|---------|
| **Runtime (CFD)** | Integrated zstd | ~10% | Automatic | Build with libzstd (recommended) |
| **Archive/Transfer** | 7z lzma2 -mx9 | ~6% | 10-11 MB/s | `7z a -m0=lzma2 -mx=9 archive.7z *.dvz` |
| **Real-time streaming** | lz4 | ~14-18% | 600-800 MB/s | `lz4 file.dvz` |
| **Quick compression** | zstd -1 | ~10-14% | 500-600 MB/s | `zstd -1 file.dvz` |

### Key Findings

1. **Best compression ratio**: 7z LZMA2 and xz achieve ~6% (94% reduction)
2. **Best speed/ratio balance**: zstd -3 at 10% ratio with 300-400 MB/s throughput
3. **zstd -3 is 40× faster than xz** with only 60% more space
4. **DVZT compresses slightly better than DVZ** due to temporal correlation patterns
5. **zstd --ultra -22 provides no benefit** over zstd -19 for this data type
6. **bzip2 and PPMd perform poorly** for CFD boundary data

### Archival Example

```bash
# Archive all boundary data with best compression
cd subset-case/constant/boundaryData/oldInternalFaces
7z a -m0=lzma2 -mx=9 ../boundaryData.7z */*.dvz */*.dvz.zstd */*.dvzt */*.dvzt.zstd

# Or with tar for faster archival (already zstd-compressed if built with libzstd)
tar -cf ../boundaryData.tar */*.dvz.zstd
```

## Limitations

- No temporal extrapolation - boundary data must cover full reconstruction time range
- Steady-state solvers not supported (requires transient simulation)
- For parallel extraction, `reconstructPar` must be run before `spaceTimeWindowInitCase`

## References

**Acknowledging this work**

If this software has been useful in your research, please consider citing:

Anton, A.-A. (2011). *"Space-Time Window Reconstruction in Parallel High Performance Numeric Simulations. Application for CFD"*, PhD Thesis, Politehnica University of Timisoara. Available at: https://dspace.upt.ro/jspui/handle/123456789/643

```bibtex
@phdthesis{Anton2011,
  author    = "Alin-Adrian Anton",
  title     = "Space-Time Window Reconstruction in High-Performance Numeric Simulations: Application for {CFD}",
  school    = "Universitatea Politehnica Timișoara",
  year      = "2011",
  month     = "November",
  type      = "PhD thesis",
  address   = "Timișoara, Romania",
  publisher = "Editura Politehnica",
  isbn      = "978-606-554-390-4",
  url       = "https://dspace.upt.ro/jspui/handle/123456789/643"
}
```

Or acknowledge: "This work made use of openfoam-spaceTimeWindow https://dev.cs.upt.ro/alin.anton/openfoam-spaceTimeWindow|https://github.com/cherubimro/openfoam-spaceTimeWindow"

**Additional references:**
- Wu, X., Zaki, T.A., Meneveau, C. (2020). *"High-fidelity inflow conditions for internal flow calculations using the boundary-layer/Reynolds-averaged Navier-Stokes approach"*, Phys. Rev. Fluids 5, 064607
- Lyn, D.A., Einav, S., Rodi, W., Park, J.-H. (1995). *"A laser-Doppler velocimetry study of ensemble-averaged characteristics of the turbulent near wake of a square cylinder"*, J. Fluid Mech. 304, 285-319
- Barry, P.J., Goldman, R.N. (1988). *"A recursive evaluation algorithm for a class of Catmull-Rom splines"*, ACM SIGGRAPH Computer Graphics 22(4), 199-204

# spaceTimeWindow Library

OpenFOAM library for space-time window extraction and reconstruction of LES simulations.

## Overview

This library provides three components for extracting a spatial subset from a full LES simulation and reconstructing the flow within that subset using pre-computed boundary conditions:

1. **spaceTimeWindowExtract** - Function object to extract boundary data during the original simulation
2. **spaceTimeWindowInitCase** - Utility to initialize a reconstruction case from extracted data
3. **spaceTimeWindow** - Boundary condition to apply extracted data during reconstruction

## Workflow

The workflow is strictly sequential:

1. **Extraction phase**: Run the original simulation with `spaceTimeWindowExtract` function object
2. **Case setup**: Run `spaceTimeWindowInitCase` to configure the reconstruction case
3. **Reconstruction phase**: Run the solver directly - everything is pre-configured

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
        fields          (U p nut);               // Fields to extract

        writeControl    timeStep;
        writeInterval   1;
    }
}
```

**Parameters:**

| Parameter    | Description                                    | Required |
|--------------|------------------------------------------------|----------|
| box          | Bounding box as ((minX minY minZ) (maxX maxY maxZ)) | yes |
| outputDir    | Output case directory for extracted data       | yes      |
| fields       | List of fields to extract                      | yes      |

**Field Selection:**

The `fields` parameter explicitly lists which fields to extract. Fields are **not** auto-detected from the turbulence model. Common configurations:

| Turbulence Model | Recommended fields |
|------------------|-------------------|
| LES Smagorinsky  | `(U p nut)` |
| LES dynamicKEqn  | `(U p nut k)` |
| RANS k-epsilon   | `(U p nut k epsilon)` |
| RANS k-omega SST | `(U p nut k omega)` |

If a requested field doesn't exist in the simulation, extraction will fail with a "field not found" error.

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
                    U               # Face-interpolated velocity
                    p               # Face-interpolated pressure
                    nut             # Face-interpolated turbulence viscosity
    <startTime>/                # Initial subset fields
        U
        p
        nut
```

**Notes:**
- Box must be **fully internal** to the domain (no intersection with external boundaries)
- The mesh is written with only `oldInternalFaces` patch (type `patch`)
- Boundary data is written at every timestep (in `execute()`), not just at write intervals
- The `extractionMetadata` file includes the list of all extracted timesteps (exact directory names)
- Currently only supports serial execution (will error if run in parallel)
- Face values are computed using linear interpolation: `U_face = w * U_inside + (1-w) * U_outside`

### spaceTimeWindowInitCase (Utility)

Initializes a fully configured reconstruction case from the extracted data.

```bash
# Run from the extraction output directory
spaceTimeWindowInitCase -sourceCase ../ufr2-02

# Or specify the extract directory explicitly
spaceTimeWindowInitCase -sourceCase ../ufr2-02 -extractDir ./subset-case
```

**Options:**

| Option         | Description                                    | Required |
|----------------|------------------------------------------------|----------|
| -sourceCase    | Source case directory (where extraction ran)   | yes      |
| -extractDir    | Directory with extracted data (default: cwd)   | no       |
| -overwrite     | Overwrite existing files                       | no       |

**What it creates:**

1. `system/controlDict` - With matching solver, deltaT, adjustTimeStep from extraction
   - **Note:** `endTime` is set to the **second-to-last** available timestep, because the `spaceTimeWindow` BC requires data from the "next" timestep for temporal interpolation
2. `system/fvSchemes`, `system/fvSolution` - Copied from source case
3. `constant/` files - All physics properties copied (mandatory for fidelity):
   - `turbulenceProperties`, `transportProperties`
   - `momentumTransport`, `thermophysicalProperties`
   - `LESProperties`, `RASProperties`, `g`
4. Initial field files with `spaceTimeWindow` BC on `oldInternalFaces`
5. Turbulence fields (nut, k, epsilon, omega, etc.) with appropriate wall functions

**After running:**
```bash
cd subset-case
pimpleFoam    # Run directly - everything is configured
```

### spaceTimeWindow (Boundary Condition)

Reads pre-computed face values from boundaryData and applies them as boundary conditions.

```cpp
oldInternalFaces
{
    type            spaceTimeWindow;
    dataDir         "constant/boundaryData";    // Path to boundaryData
    fixesValue      true;                       // See mass conservation section
    value           uniform (0 0 0);
}
```

**Properties:**

| Property       | Description                                    | Required | Default              |
|----------------|------------------------------------------------|----------|----------------------|
| dataDir        | Path to boundaryData directory                 | no       | constant/boundaryData |
| fieldTableName | Name of field file (if different from field)   | no       | field name           |
| setAverage     | Adjust mapped field to match average           | no       | false                |
| offset         | Offset value added to field                    | no       | Zero                 |
| fixesValue     | Report to adjustPhi that values are fixed      | no       | true                 |

**Notes:**
- Does NO spatial interpolation (values are pre-computed for exact face positions)
- Performs linear temporal interpolation between available timesteps if current time falls between two samples
- Will error (not extrapolate) if simulation time is outside available data range
- Metadata validation enforces identical `deltaT` to minimize interpolation effects

## Mass Conservation and `fixesValue`

The `fixesValue` option controls how `adjustPhi()` handles flux correction on the `oldInternalFaces` patch.

### fixesValue = true (default)

- `adjustPhi()` excludes this patch from flux correction
- Preserves exact boundary values from original simulation
- Any mass imbalance is corrected on other adjustable patches (typically outlets)
- Use when: You want to preserve the exact extracted boundary values

### fixesValue = false

- `adjustPhi()` includes this patch in flux correction
- Allows small modifications to the flux to ensure mass conservation
- The correction is distributed across all faces of this patch
- Use when: You want `adjustPhi()` to heal/repair mass imbalance caused by face interpolation

**Recommendation:** For LES reconstruction where accurate inflow physics matter, `fixesValue = false` may be preferable as it distributes the mass correction uniformly across the `oldInternalFaces` rather than dumping it all on the outlet.

## Time Settings Validation

The library enforces identical time settings between extraction and reconstruction to prevent temporal interpolation/extrapolation errors.

### Metadata Storage

During extraction, `spaceTimeWindowExtract` writes an `extractionMetadata` file containing:

```
constant/boundaryData/oldInternalFaces/extractionMetadata
```

Contents:
- `openfoamVersion` - OpenFOAM version string (e.g., "v2512")
- `openfoamApi` - OpenFOAM API number (e.g., 2512)
- `solver` - Application/solver name used during extraction (e.g., "pimpleFoam")
- `deltaT` - Time step used during extraction
- `adjustTimeStep` - Whether adaptive time stepping was enabled
- `timePrecision` - Floating point precision for time values
- `extractionStartTime` - Start time of extraction
- `fixedDeltaT` - Fixed deltaT from controlDict (if specified)

### Validation During Reconstruction

The `spaceTimeWindow` BC validates these settings on first use:

1. **OpenFOAM API mismatch** - Warning if reconstruction uses different OpenFOAM version
2. **Solver mismatch** - Warning if reconstruction uses different solver
3. **deltaT mismatch** - Fatal error if reconstruction uses different time step
4. **adjustTimeStep mismatch** - Fatal error if adaptive time stepping setting differs
5. **adjustTimeStep warning** - Warning if both use adaptive time stepping (time steps may not align exactly)

Example error message:
```
FOAM FATAL ERROR:
Time step mismatch between extraction and reconstruction!
    Extraction deltaT: 8.53e-04
    Current deltaT:    1e-03
    Patch: oldInternalFaces

    The spaceTimeWindow BC requires identical time settings
    to avoid temporal interpolation/extrapolation errors.
    Please set deltaT = 8.53e-04 in system/controlDict
```

**Best Practice:** Use fixed time stepping (`adjustTimeStep = no`) with identical `deltaT` values for both extraction and reconstruction.

## Building

```bash
# Build the library (BC and function object)
cd src/spaceTimeWindow
wmake

# Build the utility
cd applications/utilities/preProcessing/spaceTimeWindowInitCase
wmake
```

## Requirements

- OpenFOAM v2512 or compatible version
- The subset mesh must be created using `subsetMesh` with exact cell extraction (no interpolation)

## RANS Compatibility

The library works with transient RANS simulations (e.g., `pimpleFoam` with k-epsilon or k-omega SST):

**Supported:**
- Transient RANS with time-varying boundary conditions
- All transported turbulence variables (k, epsilon, omega)
- Wall function handling (though subset typically has no walls)

**Not supported:**
- Steady-state RANS (`simpleFoam`) - the "space-time" concept requires time evolution

**Example for k-omega SST:**
```cpp
fields    (U p nut k omega);
```

The extraction and reconstruction process is identical to LES. Ensure all transported turbulence quantities are included in the `fields` list.

## Limitations

- Parallel execution not supported for extraction (run on reconstructed case)
- No spatial interpolation - mesh topology must match exactly
- Temporal extrapolation not allowed - boundary data must cover full simulation time
- Steady-state solvers not supported (requires transient simulation)

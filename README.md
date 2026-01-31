# spaceTimeWindow Library

A library for space-time window extraction and reconstruction of LES simulations, developed for OpenFOAM.

**Tested with:** OpenFOAM v2512 (openfoam.com)

> **Trademark Notice:** This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM and OpenCFD trade marks. OPENFOAM is a registered trade mark of OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com.

## Overview

This library provides three components for extracting a spatial subset from a full LES simulation and reconstructing the flow within that subset using pre-computed boundary conditions:

1. **spaceTimeWindowExtract** - Function object to extract boundary data during the original simulation
2. **spaceTimeWindowInitCase** - Utility to initialize a reconstruction case from extracted data
3. **spaceTimeWindow** - Boundary condition to apply extracted data during reconstruction

## Workflow

The workflow is strictly sequential:

### Serial Execution

1. **Extraction phase**: Run the original simulation with `spaceTimeWindowExtract` function object
2. **Case setup**: Run `spaceTimeWindowInitCase` to configure the reconstruction case
3. **Reconstruction phase**: Run the solver directly - everything is pre-configured

### Parallel Execution

1. **Extraction phase**: Run the original simulation in parallel with `spaceTimeWindowExtract` function object
   - The extractor automatically forces a field write at the extraction start time
   - Boundary data is gathered from all processors and written by master
2. **Reconstruct fields**: Run `reconstructPar -time <extractionStartTime>` to reconstruct the initial time
3. **Case setup**: Run `spaceTimeWindowInitCase -sourceCase <path>` to create subset mesh and initial fields
4. **Reconstruction phase**: Run the solver on the subset case

```bash
# Example parallel workflow
cd source-case
mpirun -np 4 pimpleFoam -parallel    # Extraction happens automatically

# Reconstruct the extraction start time
reconstructPar -time 0.0001          # Use actual extraction start time

# Initialize the subset case
cd ../subset-case
spaceTimeWindowInitCase -sourceCase ../source-case

# Run reconstruction (can be serial or parallel)
pimpleFoam
```

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

        // Write format for boundary data (optional)
        // Options: ascii, binary, deltaVarint
        // Default: ascii
        writeFormat     deltaVarint;

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
| fields       | List of fields to extract                      | yes      | - |
| writeFormat  | Format for boundary data: `ascii`, `binary`, or `deltaVarint` | no | ascii |
| deltaVarintPrecision | Decimal digits for delta-varint quantization | no | 6 |
| writeCompression | Gzip compression for ascii/binary (ignored for deltaVarint) | no | off |

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
                    U               # Face-interpolated velocity (or U.dvz)
                    p               # Face-interpolated pressure (or p.dvz)
                    nut             # Face-interpolated turbulence viscosity (or nut.dvz)
    <startTime>/                # Initial subset fields (always ASCII)
        U
        p
        nut
```

**Note:** With `writeFormat deltaVarint`, boundary data files have `.dvz` extension (e.g., `U.dvz`). Initial fields remain ASCII regardless of writeFormat setting.

**Notes:**
- Box must be **fully internal** to the domain (no intersection with external boundaries)
- The mesh is written with only `oldInternalFaces` patch (type `patch`)
- Boundary data is written at every timestep (in `execute()`), not just at write intervals
- The `extractionMetadata` file includes the list of all extracted timesteps (exact directory names)
- **Supports both serial and parallel execution**
- Face values are computed using linear interpolation: `U_face = w * U_inside + (1-w) * U_outside`

**Parallel Execution:**
- The extraction box can span multiple processor domains
- Boundary data from all processors is gathered to the master and written as a single file
- Processor boundary faces (where the extraction boundary crosses processor boundaries) are handled automatically
- Field interpolation correctly exchanges values between processors for faces on processor boundaries
- **Automatic field write**: At extraction start time, the extractor forces the solver to write all fields to disk
- **Mesh handling**: In parallel mode, only extraction box parameters are written (not the mesh)
- `spaceTimeWindowInitCase` creates the subset mesh from the reconstructed source case, ensuring correct cell ordering

**Parallel Extraction Output:**
```
outputDir/
    constant/
        polyMesh/
            extractionBox       # Extraction parameters for mesh creation
        boundaryData/
            oldInternalFaces/
                points              # Face centres (gathered from all processors)
                extractionMetadata  # Settings, timestep list, nProcs info
                <time>/
                    U               # Gathered face-interpolated fields (or U.dvz)
                    p               # (or p.dvz with deltaVarint)
                    nut             # (or nut.dvz)
```

**Important:** After parallel extraction, you must run `reconstructPar -time <startTime>` before `spaceTimeWindowInitCase`.

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
   - **Note:** `startTime` is set to the **third** available timestep (t_2), because cubic interpolation needs t_0 and t_1 as lookback data
   - **Note:** `endTime` is set to the **third-to-last** available timestep (t_{n-2}), because cubic interpolation needs t_{n-1} and t_n as lookahead data
2. `system/fvSchemes`, `system/fvSolution` - Copied from source case
3. `constant/` files - All physics properties copied (mandatory for fidelity):
   - `turbulenceProperties`, `transportProperties`
   - `momentumTransport`, `thermophysicalProperties`
   - `LESProperties`, `RASProperties`, `g`
4. Initial field files with `spaceTimeWindow` BC on `oldInternalFaces`
5. Turbulence fields (nut, k, epsilon, omega, etc.) with appropriate wall functions

**Parallel case handling:**

When the extraction was done in parallel (detected by `extractionBox` file):
- Creates subset mesh from the reconstructed source case using the extraction bounding box
- Extracts initial fields from the source case at the extraction start time
- Requires that `reconstructPar -time <startTime>` was run on the source case first

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
| allowTimeInterpolation | Permit interpolation for missing timesteps | no | false |
| timeInterpolationScheme | Interpolation method: `linear` or `cubic` | no | linear |
| reportFlux     | Print net flux through patch (velocity only)   | no       | false                |

**Notes:**
- Does NO spatial interpolation (values are pre-computed for exact face positions)
- By default (`allowTimeInterpolation = false`), requires exact timestep matching with source simulation
- Will error if simulation time doesn't match an available timestep (unless `allowTimeInterpolation = true`)
- Metadata validation enforces identical `deltaT` to ensure timestep alignment
- **Extrapolation outside the extraction window is NEVER allowed**, regardless of `allowTimeInterpolation` setting

### Time Interpolation

By default, the BC requires exact timestep matching between extraction and reconstruction. This is the safest option ensuring reproducibility.

If your reconstruction case cannot match the extraction timesteps exactly (e.g., different solver constraints, restarts), you can enable time interpolation:

```cpp
oldInternalFaces
{
    type                    spaceTimeWindow;
    dataDir                 "constant/boundaryData";
    allowTimeInterpolation  true;   // Permit interpolation between timesteps
    timeInterpolationScheme cubic;  // Use cubic spline (or "linear")
    value                   uniform (0 0 0);
}
```

**When `allowTimeInterpolation = false` (default):**
- Simulation time must exactly match an available sample time (within 1% of deltaT tolerance)
- Fatal error if no matching timestep is found
- Ensures bit-reproducible results

**When `allowTimeInterpolation = true`:**
- Interpolation between bracketing timesteps is permitted
- Useful when extraction and reconstruction have slightly different time stepping
- Extrapolation outside `[extractionStartTime, extractionStopTime]` still causes a fatal error

**Time Interpolation Schemes:**

| Scheme | Description | Points Used |
|--------|-------------|-------------|
| `linear` | Simple linear interpolation (default) | 2 (bracketing timesteps) |
| `cubic` | Catmull-Rom cubic spline | 4 (2 before, 2 after) |

**Linear interpolation** (`timeInterpolationScheme linear`):
- Uses two bracketing timesteps: t_i and t_{i+1}
- Simple weighted average: `value = (1-alpha) * v_i + alpha * v_{i+1}`
- Sufficient for smoothly varying flows

**Cubic spline interpolation** (`timeInterpolationScheme cubic`):
- Uses four timesteps: t_{i-1}, t_i, t_{i+1}, t_{i+2}
- **Centripetal Catmull-Rom spline** - handles non-uniform time spacing correctly
- Essential for adaptive timestepping (`adjustTimeStep = yes`) where dt varies with CFL
- Provides C1 continuity (smooth first derivatives) without overshoots or cusps
- Better captures temporal variations in turbulent flows
- Requires at least 5 timesteps in boundaryData (2 buffer at start + 2 buffer at end + 1 usable)

The centripetal parameterization (Barry & Goldman, 1988) prevents the numerical artifacts that standard Catmull-Rom produces with non-uniform data spacing, making it the correct choice for CFD simulations with variable timesteps.

**Time range for cubic interpolation:** `spaceTimeWindowInitCase` automatically sets:
- `startTime = t_2` (third timestep) - so t_0 and t_1 are available for lookback
- `endTime = t_{n-2}` (third-to-last) - so t_{n-1} and t_n are available for lookahead

This ensures cubic interpolation has all 4 required points throughout the entire reconstruction. When the solver is at time t in the interval [t_2, t_3], cubic interpolation uses t_1, t_2, t_3, t_4 - meaning t_0 and t_1 provide the necessary lookback buffer. Requires at least 5 timesteps in boundaryData.

For turbulent LES simulations with significant temporal fluctuations, `cubic` interpolation produces smoother, more physically realistic boundary conditions when interpolation is needed

### Flux Reporting

The `reportFlux` option enables diagnostic output showing the net mass flux through the `oldInternalFaces` patch at each timestep. This is useful for verifying mass conservation.

```cpp
oldInternalFaces
{
    type            spaceTimeWindow;
    dataDir         "constant/boundaryData";
    fixesValue      true;
    reportFlux      true;   // Enable flux diagnostics
    value           uniform (0 0 0);
}
```

**Output example:**
```
spaceTimeWindow flux report [oldInternalFaces] t=0.001 fixesValue=true netFlux=1.234e-06 (in=-0.0523 out=0.0523)
```

- **netFlux**: Total net flux through patch (should be ~0 for mass conservation)
- **in**: Sum of inward fluxes (negative values, into domain)
- **out**: Sum of outward fluxes (positive values, out of domain)
- **fixesValue**: Shows whether `adjustPhi()` will modify this patch's flux

When `fixesValue=true`, `adjustPhi()` treats the patch as fixed and won't correct its flux - any imbalance is absorbed by other patches (typically outlets). When `fixesValue=false`, `adjustPhi()` can adjust the flux on this patch to enforce global mass conservation.

**Note:** `reportFlux` only produces output for vector fields (i.e., velocity U).

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

## Boundary Data Compression

The `writeFormat` parameter controls how boundary data files are written. For long LES simulations with many timesteps, storage can become significant.

### Format Comparison

| Format | Extension | Typical Size | Notes |
|--------|-----------|--------------|-------|
| `ascii` | (none) | 100% baseline | Human-readable, default |
| `binary` | (none) | ~50% | OpenFOAM native binary |
| `ascii` + gzip | .gz | ~10% | `writeCompression on` |
| `binary` + gzip | .gz | ~8% | `writeFormat binary` + `writeCompression on` |
| `deltaVarint` | .dvz | **~2.7%** | Best compression, recommended |

### ASCII and Binary Formats

Standard OpenFOAM formats with optional gzip compression:

```cpp
// ASCII (default, human-readable)
writeFormat     ascii;

// Binary (smaller, faster I/O)
writeFormat     binary;

// Either format with gzip compression
writeCompression on;
```

### Delta-Varint Codec (`writeFormat deltaVarint`)

A specialized codec optimized for time-series CFD data that achieves ~97% storage reduction:

**How it works:**
1. **Component-major ordering**: Stores all Ux values, then Uy, then Uz (groups similar values)
2. **Delta encoding**: Stores differences between consecutive values (small deltas)
3. **Quantization**: Rounds deltas to configurable precision (default: 6 decimal digits)
4. **Varint encoding**: Variable-length integer encoding (1-9 bytes based on magnitude)
5. **Zigzag encoding**: Efficient encoding of signed integers

**File format:**
- Magic number: `DVZ1` (0x315A5644)
- Header: element count, component count, precision
- First value: stored as raw 8-byte double (exact)
- Subsequent values: delta-encoded varints

**Precision:**
- `deltaVarintPrecision 6` means deltas quantized to `round(delta * 10^6)`
- Precision loss is ~1e-6 relative to consecutive value differences
- For typical CFD data, this is well below solver tolerance

**What uses deltaVarint:**
- Boundary data files only (`constant/boundaryData/.../U.dvz`, `p.dvz`)

**Important:** When using `writeFormat deltaVarint`, all other output is written as **uncompressed ASCII**:
- Initial fields (`<startTime>/U`, `p`) - ASCII
- Subset mesh (`constant/polyMesh/`) - ASCII
- Metadata (`extractionMetadata`) - ASCII

The `writeCompression` setting is **ignored** when using deltaVarint. This is intentional because:
1. Boundary data dominates storage (written every timestep)
2. Initial fields and mesh are written once
3. The ~2.7% compression ratio already far exceeds gzip

**Auto-detection:**
The `spaceTimeWindow` BC and `spaceTimeWindowInitCase` automatically detect `.dvz` files and read them correctly. No configuration needed on the reconstruction side.

## Encryption (Optional)

Boundary data can be encrypted using X25519 asymmetric encryption (libsodium sealed boxes). This allows extraction with a public key while requiring the private key for reconstruction.

### Building with Encryption Support

```bash
cd src/spaceTimeWindow
export FOAM_USE_SODIUM=1 && wmake

cd applications/utilities/preProcessing/spaceTimeWindowInitCase
export FOAM_USE_SODIUM=1 && wmake

cd ../spaceTimeWindowKeygen
export FOAM_USE_SODIUM=1 && wmake
```

Requires libsodium development package (`libsodium-dev` on Debian/Ubuntu, `libsodium-devel` on RHEL/Fedora).

### Key Generation

```bash
spaceTimeWindowKeygen
# Output:
# Public key:  fqzYQ0U8j27tFEr5WzEMylbvXYP+9CAyk0JhwwZ2rwg=
# Private key: QgzxB5b+DGPQH8exbWDe18n4Kv0nu5gqljI2RPBCwl4=
#
# IMPORTANT: Store the private key securely!
```

### Encrypted Extraction

Add the public key to the function object configuration:

```cpp
functions
{
    extractSubset
    {
        type            spaceTimeWindowExtract;
        libs            (spaceTimeWindow);

        box             ((0.05 -0.25 0.01) (0.90 0.25 0.38));
        outputDir       "../subset-case";
        fields          (U p nut);

        writeFormat     deltaVarint;
        deltaVarintPrecision  7;

        // Encryption with X25519 public key
        publicKey       "fqzYQ0U8j27tFEr5WzEMylbvXYP+9CAyk0JhwwZ2rwg=";

        writeControl    timeStep;
        writeInterval   1;
    }
}
```

Encrypted files have `.enc` extension (e.g., `U.dvz.enc` for compressed+encrypted).

### Decryption During Case Initialization

```bash
cd subset-case
spaceTimeWindowInitCase -sourceCase ../source-case
# Prompts for private key (no echo):
# Enter private key (base64):
```

The public key is automatically derived from the private key. All encrypted boundary data is decrypted in-place during initialization.

### Security Properties

- **Sealed box encryption**: Anonymous sender, only recipient (private key holder) can decrypt
- **X25519 key derivation**: Public key can be computed from private key
- **No unencrypted data on disk**: Encryption happens before writing during extraction
- **Compile-time optional**: Without `FOAM_USE_SODIUM=1`, encryption code is not included

## Time Settings Validation

The library enforces identical time settings between extraction and reconstruction to ensure exact timestep matching.

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
    for exact timestep matching.
    Please set deltaT = 8.53e-04 in system/controlDict
```

**Best Practice:** Use fixed time stepping (`adjustTimeStep = no`) with identical `deltaT` values for both extraction and reconstruction. All timesteps from the source simulation are required.

## Building

```bash
# Build everything (library + utilities)
./Allwmake

# Or build with encryption support (requires libsodium)
# ./Allwmake will auto-detect libsodium if installed
```

### Manual Build

```bash
# Build the library (BC and function object)
cd src/spaceTimeWindow
wmake

# Build the utility
cd applications/utilities/preProcessing/spaceTimeWindowInitCase
wmake
```

## Documentation

API documentation can be generated with Doxygen:

```bash
./Allwmake doc
```

This generates HTML documentation in `doc/html/index.html`.

**Requirements:** Doxygen and Graphviz (for class diagrams)

## Requirements

- **OpenFOAM v2512** (openfoam.com) or compatible ESI-OpenCFD version
- The subset mesh must be created using `subsetMesh` with exact cell extraction (no interpolation)
- For parallel extraction: MPI environment for distributed execution

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

## Example Case

The `examples/ufr2-02` directory contains a complete LES test case:

- **Case**: ERCOFTAC UFR2-02 - Flow around a square cylinder at Re = 22,000
- **Acknowledgment**: Mesh generation script by Niklas Nordin (ERCOFTAC Classic Collection Database Case 043)
- **Reference**: Lyn et al. (1995) J. Fluid Mech. 304, 285-319

```bash
cd examples/ufr2-02
./Allrun
```

See `examples/ufr2-02/README.md` for detailed instructions.

## References

**Acknowledging this work**

If this software has been useful in your research, please consider:
- Citing: 
Anton, A.-A. (2011). *"Space-Time Window Reconstruction in Parallel High Performance Numeric Simulations. Application for CFD"*, PhD Thesis, Politehnica University of Timisoara. Available at: https://dspace.upt.ro/jspui/handle/123456789/643

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

- Or acknowledging: "This work made use of openfoam-spaceTimeWindow https://dev.cs.upt.ro/alin.anton/openfoam-spaceTimeWindow"

**Additional references:**
- Lyn, D.A., Einav, S., Rodi, W., Park, J.-H. (1995). *"A laser-Doppler velocimetry study of ensemble-averaged characteristics of the turbulent near wake of a square cylinder"*, J. Fluid Mech. 304, 285-319
- ERCOFTAC Classic Collection Database: Case 043

## Limitations

- No spatial interpolation - mesh topology must match exactly
- No temporal extrapolation - boundary data must cover full reconstruction time range
- By default, no temporal interpolation (requires exact timestep matching); can be enabled via `allowTimeInterpolation` with `linear` or `cubic` scheme
- Steady-state solvers not supported (requires transient simulation)
- For parallel extraction, `reconstructPar` must be run before `spaceTimeWindowInitCase`

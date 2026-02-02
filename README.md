# spaceTimeWindow Library

A library for space-time window extraction and reconstruction of LES simulations, developed for OpenFOAM.

**Tested with:** OpenFOAM v2512 (openfoam.com)

> **Trademark Notice:** This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM and OpenCFD trade marks. OPENFOAM is a registered trade mark of OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com.

## Overview

This library provides three components for extracting a spatial subset from a full LES simulation and reconstructing the flow within that subset using pre-computed boundary conditions:

1. **spaceTimeWindowExtract** - Function object to extract boundary data during the original simulation
2. **spaceTimeWindowInitCase** - Utility to initialize a reconstruction case from extracted data
3. **spaceTimeWindow** / **spaceTimeWindowInletOutlet** - Boundary conditions to apply extracted data during reconstruction

## Key Concepts

### Boundary Condition Approaches

The library provides two approaches for applying boundary conditions on the extraction boundary (`oldInternalFaces`):

#### 1. Inlet-Outlet BC (`-inletOutletBC`) - **Recommended for Unsteady Flows**

Uses `spaceTimeWindowInletOutlet` BC which applies:
- **Velocity (U)**: Flux-based switching - Dirichlet (prescribed value) at inflow faces, zeroGradient at outflow faces
- **All scalar fields (p, nut, k, epsilon, omega)**: zeroGradient

This approach is **physically correct** for unsteady turbulent flows because:
- Turbulent flows have instantaneous velocity fluctuations that can reverse direction locally
- Vortex shedding, recirculation zones, and turbulent eddies cause portions of the boundary to alternate between inflow and outflow
- Prescribing fixed values at outflow faces is non-physical (information should leave the domain freely)
- The flux-based switching automatically adapts to the instantaneous flow direction at each face

```bash
spaceTimeWindowInitCase -sourceCase ../source-case -inletOutletBC
```

#### 2. Fixed Outlet Direction (`-outletDirection`) - For Steady-Mean Flows

Creates a separate `outlet` patch from faces at one edge of the extraction box:
- **oldInternalFaces**: Pure Dirichlet BC (spaceTimeWindow) on all remaining faces
- **outlet**: Pressure relief patch with inletOutlet for U and zeroGradient for p

Use when the mean flow direction is well-defined and outflow always occurs at a known location.

```bash
spaceTimeWindowInitCase -sourceCase ../source-case -outletDirection "(1 0 0)"
```

**Note:** These options are mutually exclusive. Using both together produces an error with guidance.

### Initial Fields vs Boundary Data Fields

The extraction can separate which fields are used for:
- **Initial conditions** (`initialFields`): Fields written to the start time directory for solver initialization
- **Time-varying boundary data** (`fields`): Fields written to `boundaryData/` for time-varying BCs

This allows extracting more fields for initial conditions (e.g., `U p nut k omega`) while only storing time-varying data for velocity:

```cpp
// In extraction function object
fields          (U);           // Only U for time-varying BC (saves storage)
initialFields   (U p nut);     // More fields for initial conditions
```

Fields NOT in `boundaryData` automatically get `zeroGradient` BC on `oldInternalFaces`.

## Workflow

### Serial Execution

```bash
# 1. Run extraction during simulation
cd source-case
pimpleFoam    # With spaceTimeWindowExtract function object

# 2. Initialize reconstruction case (recommended: inlet-outlet BC)
cd ../subset-case
spaceTimeWindowInitCase -sourceCase ../source-case -inletOutletBC

# 3. Run reconstruction
pimpleFoam
```

In serial mode:
- Mesh and initial fields are written directly to the output directory
- Boundary data is written at every timestep
- The extractor forces field writes at t_1 and t_2 (for interpolation buffers)

### Parallel Execution

```bash
# 1. Run extraction during parallel simulation
cd source-case
mpirun -np 4 pimpleFoam -parallel    # Extraction happens automatically

# 2. Reconstruct the start time (t_2 for cubic interpolation)
reconstructPar -time 0.0002

# 3. Initialize reconstruction case (recommended: inlet-outlet BC)
cd ../subset-case
spaceTimeWindowInitCase -sourceCase ../source-case -inletOutletBC

# 4. Run reconstruction (serial or parallel)
pimpleFoam
```

In parallel mode:
- Boundary data from all processors is gathered to master and written as single files
- The extractor forces field writes at t_0, t_1, and t_2 for interpolation buffers (t_0 only in parallel)
- Only extraction box parameters are written (not mesh) - `spaceTimeWindowInitCase` creates the mesh
- `reconstructPar` must be run at the reconstruction start time (t_2) to provide source fields

### Time Interpolation Modes

| Mode | Start Time | Timesteps Used | Use Case |
|------|------------|----------------|----------|
| `none` (exact) | t_0 | 1 (exact match) | "Bit"-reproducible results |
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
        fields          (U);                     // Typically just velocity

        // Fields for initial conditions (optional, defaults to 'fields')
        initialFields   (U p nut);               // More fields for IC

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

For unsteady turbulent flows with `-inletOutletBC`, the recommended approach is:

```cpp
fields          (U);               // Only velocity needs time-varying BC data
initialFields   (U p nut);         // All fields needed for solver startup
```

This works because:
- Velocity is the only field that needs prescribed values at inflow faces
- Pressure and turbulence quantities use zeroGradient (no BC data needed)
- Reduces storage by ~66% compared to extracting all fields

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
    <startTime>/                # Initial subset fields (always ASCII)
        U
        p
        nut
```

### spaceTimeWindowInitCase (Utility)

Initializes a fully configured reconstruction case from the extracted data.

```bash
# Recommended: inlet-outlet BC for unsteady turbulent flows
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
| -inletOutletBC   | **Recommended.** Use flux-based inlet-outlet BC for U, zeroGradient for scalars | no |
| -outletDirection | Create fixed outlet patch in given direction (e.g., "(1 0 0)") | no |
| -outletFraction  | Fraction of box extent for outlet region (default: 0.1) | no |
| -correctMassFlux | Apply least-squares mass flux correction to boundaryData | no |
| -initialFields   | Override initial fields list (e.g., "(U p nut k)") | no |
| -overwrite       | Overwrite existing files                         | no       |

**Note:** `-inletOutletBC` and `-outletDirection` are mutually exclusive.

**What it creates:**

1. `system/controlDict` - With matching solver, deltaT, adjustTimeStep from extraction
   - `startTime` set to t_2 (third timestep) for cubic interpolation buffer
   - `endTime` set to t_{n-2} (third-to-last) for cubic interpolation buffer
2. `system/fvSchemes`, `system/fvSolution` - Copied from source case (with pRefPoint added)
3. `constant/` files - All physics properties copied (mandatory for fidelity)
4. Initial field files with appropriate BCs based on options and boundaryData availability

**Boundary Condition Assignment Logic:**

For each field in the initial time directory:

| Field in boundaryData? | -inletOutletBC? | Field Type | BC Applied |
|------------------------|-----------------|------------|------------|
| Yes | Yes | vector (U) | `spaceTimeWindowInletOutlet` |
| Yes | Yes | scalar | `zeroGradient` |
| Yes | No | any | `spaceTimeWindow` (Dirichlet) |
| No | any | any | `zeroGradient` |

This ensures fields without time-varying data automatically get appropriate BCs.

### spaceTimeWindowInletOutlet (Boundary Condition) - **Recommended**

Flux-based boundary condition that reads pre-computed velocity values and applies them only at inflow faces.

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

**1. Use `-inletOutletBC` (Recommended)**

The flux-based inlet-outlet BC naturally handles mass conservation:
- Outflow faces use zeroGradient, allowing natural outflow
- No artificial mass imbalance from prescribed outflow velocities
- Works without any special mass correction

**2. Use `-correctMassFlux`**

Applies least-squares correction to boundaryData to ensure exact mass conservation:

```bash
spaceTimeWindowInitCase -sourceCase ../source -inletOutletBC -correctMassFlux
```

The correction minimizes ||U_corrected - U||² subject to Σ(U_corrected · Sf) = 0:
```
U_corrected = U - (imbalance / totalSfMag) * n
```

**3. Use `-outletDirection` with `fixesValue true`**

Creates an outlet patch where mass imbalance can escape:
- `oldInternalFaces` uses `fixesValue true` (values not modified by adjustPhi)
- `outlet` uses `inletOutlet` BC (allows adjustPhi correction)

### Recommended Configurations

| Use Case | Command | Notes |
|----------|---------|-------|
| Unsteady turbulent (vortex shedding, etc.) | `-inletOutletBC` | **Recommended**. Natural mass balance. |
| Unsteady turbulent + extra safety | `-inletOutletBC -correctMassFlux` | Ensures exact mass balance. |
| Steady-mean flow direction | `-outletDirection "(1 0 0)"` | When outlet location is known. |
| Maximum fidelity | `-correctMassFlux` (no outlet) | All faces prescribed, exact mass balance. |

## Boundary Data Compression

The `writeFormat` parameter controls how boundary data files are written.

| Format | Extension | Typical Size | Notes |
|--------|-----------|--------------|-------|
| `ascii` | (none) | 100% baseline | Human-readable, default |
| `binary` | (none) | ~50% | OpenFOAM native binary |
| `ascii` + gzip | .gz | ~10% | `writeCompression on` |
| `binary` + gzip | .gz | ~8% | Combined |
| `deltaVarint` | .dvz | **~2.7%** | High compression, self-contained |
| `dvzt` | .dvzt | **~2.4%** | Best compression, recommended |

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

During extraction, DVZT writes smaller .dvzt files. During case initialization, `spaceTimeWindowInitCase` automatically converts .dvzt files to .dvz format (required because delta frames need sequential processing). The resulting .dvz files are read by the spaceTimeWindow BC at runtime.

```
Extraction:  sim.dvzt files (~10% smaller than .dvz)
     |
     v
spaceTimeWindowInitCase:  Converts .dvzt -> .dvz (sequential)
     |
     v
Runtime:  Reads .dvz files (no changes to BC code)
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

Boundary data can be encrypted using X25519 asymmetric encryption (libsodium sealed boxes).

### Building with Encryption Support

```bash
export FOAM_USE_SODIUM=1 && ./Allwmake
```

### Usage

```bash
# Generate key pair
spaceTimeWindowKeygen
# Public key:  fqzYQ0U8j27tFEr5WzEMylbvXYP+9CAyk0JhwwZ2rwg=
# Private key: QgzxB5b+DGPQH8exbWDe18n4Kv0nu5gqljI2RPBCwl4=

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
./Allwmake

# Or with encryption support (auto-detects libsodium)
./Allwmake
```

### Manual Build

```bash
# Build the library
cd src/spaceTimeWindow
wmake

# Build utilities
cd applications/utilities/preProcessing/spaceTimeWindowInitCase
wmake
```

## Documentation

API documentation can be generated with Doxygen:

```bash
./Allwmake doc
```

## Requirements

- **OpenFOAM v2512** (openfoam.com) or compatible ESI-OpenCFD version
- For parallel extraction: MPI environment
- For encryption: libsodium development package

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

## Limitations

- No spatial interpolation - mesh topology must match exactly
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

Or acknowledge: "This work made use of openfoam-spaceTimeWindow https://dev.cs.upt.ro/alin.anton/openfoam-spaceTimeWindow"

**Additional references:**
- Lyn, D.A., Einav, S., Rodi, W., Park, J.-H. (1995). *"A laser-Doppler velocimetry study of ensemble-averaged characteristics of the turbulent near wake of a square cylinder"*, J. Fluid Mech. 304, 285-319
- Barry, P.J., Goldman, R.N. (1988). *"A recursive evaluation algorithm for a class of Catmull-Rom splines"*, ACM SIGGRAPH Computer Graphics 22(4), 199-204

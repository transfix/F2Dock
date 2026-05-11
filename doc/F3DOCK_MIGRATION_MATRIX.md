# F3Dock to F2Dock Migration Matrix

This matrix captures a second-pass triage of F3Dock capabilities for incremental porting into modern F2Dock.

## Naming and runtime mode

The "**F2**" in F2Dock stands for **Fast Fourier**: the original F2Dock pipeline performs FFT-accelerated *rigid*-body protein-protein docking. The "**F3**" in F3Dock stands for **Flexible Fast Fourier**: F3Dock layers protein-flexibility sampling stages (tripeptide loop closure, hinge / shear kinematics, domain-graph composition, ICP refinement) on top of the same FFT scoring core.

Because both pipelines share the FFT scoring infrastructure, this repository ships a **single consolidated executable** (`F2Dock`) that can run in either mode at runtime. The mode is selected with a flag rather than being baked in at build time so a single binary can serve both workflows:

| Mode value (`PARAMS_IN::dockMode`) | Enum (`f3dock::DockMode`) | Pipeline |
|---|---|---|
| `0` | `kF2Dock` | F2Dock — rigid Fast Fourier docking (default) |
| `1` | `kF3Dock` | F3Dock — Flexible Fast Fourier docking |

Selecting the mode (CLI overrides input file when both are present):

```bash
F2Dock --mode=f2 parameterFile     # rigid (synonyms: f2dock, rigid)
F2Dock --mode=f3 parameterFile     # flexible (synonyms: f3dock, flex, flexible)
F2Dock --f2dock   parameterFile    # shorthand
F2Dock --f3dock   parameterFile    # shorthand
```

In the input file:

```
dockMode f2     # or: f2dock, rigid
dockMode f3     # or: f3dock, flex, flexible
```

The active mode is printed at startup (`Docking mode: F2Dock (rigid)` / `Docking mode: F3Dock (flexible)`). Defining and parsing the value lives in [inc/f3dock/DockMode.h](inc/f3dock/DockMode.h); unit tests in [tests/unit/test_dock_mode.cpp](tests/unit/test_dock_mode.cpp).

All F3Dock-derived stages (`f3dock-loop-closure`, `f3dock-flex`, `f3dock-icp`, `f3dock-domain`) are gated to run **only** when `dockMode == kF3Dock`. Existing F2Dock workflows behave identically with no input changes.

## Capability matrix

| Capability | Scientific Value | Effort | Risk | Recommended Action |
|---|---|---:|---:|---|
| Tripeptide loop closure core (Coutsias polynomial solver) | High for local flexibility and constrained loop sampling | Medium | Medium | Port now as optional module (`f3dock-loop-closure`) with unit tests |
| Hinge bending / shear kinematics (`ProteinFlexibility`) | High for domain motion exploration | High | High | Phase 2: port behind feature flag after loop-closure API stabilizes |
| Domain graph modeling (`DomainComplex`) | Medium-High for flexible docking workflows | High | High | Phase 3 candidate; requires architecture adaptation |
| ICP alignment (`libicp`) | Medium for refinement and re-scoring | Medium | Medium | Prototype as independent utility library in Phase 2 |
| Blurmaps / GOA molecular surface stack | Medium | High | High | Defer until baseline flexible pipeline lands |
| PRGN / extra geometry stack | Medium | High | High | Defer; extract only needed pieces once use-cases are clear |
| NFFT-based acceleration paths | Medium-High (for some workloads) | High | High | Keep optional; add dependency detection first, do not hard-require |
| LAPACK-backed math extensions | Medium | Low-Medium | Low-Medium | Enable optional linking in science stack now |
| Eigen-based linear algebra utilities | Medium | Low-Medium | Low | Enable optional linking in science stack now |
| Full F3Dock executable behavior parity | Very High long-term | Very High | Very High | Multi-phase roadmap, not a single PR |

## Porting Principles

1. Keep legacy-derived code optional and isolated behind build options.
2. Add deterministic unit tests before integrating solver outputs into docking flows.
3. Maintain current F2Dock integration regression suite as a release gate.
4. Introduce new scientific dependencies as optional, never mandatory.

## Current Phase (implemented in this branch)

- Added experimental `f3dock-loop-closure` static library.
- Ported tripeptide loop-closure core solver into isolated module.
- Added stable wrapper API in `inc/f3dock/loop_closure/LoopClosureSolver.h`.
- Added unit tests for feasible and infeasible closure constraints.
- Added optional scientific stack discovery (`Eigen3`, `LAPACK`, `NFFT`).

## Phase 2 Progress (current incremental port)

- Added experimental `f3dock-icp` static library (point-to-point ICP prototype).
- Imported `libicp` core components (`matrix`, `kdtree`, point-to-point ICP loop).
- Added wrapper API in `inc/f3dock/icp/IcpAligner.h`.
- Added deterministic unit tests for rigid alignment recovery and input validation.
- Added experimental `f3dock-flex` scaffold for hinge/shear kinematics helpers.
- Added deterministic unit tests for hinge rotation and planar shear behavior.
- Added experimental `f3dock-domain` scaffold for domain-graph transform composition.
- Added deterministic unit tests for world-transform chaining and cycle rejection.

## Phase 2 Completion (libicp parity)

- Imported point-to-plane ICP (`libicp::IcpPointToPlane`) including
  PCA-based normal estimation from the model point cloud.
- Extended `IcpAligner` with `alignPointToPlane(...)` mirroring the
  point-to-point API, gated by a minimum neighbor count (>= 3).
- Added unit tests for surface-fit convergence and input validation
  (too-few points, too-small neighborhoods).

## Phase 3.5 (libcvc volume I/O adoption)

- Replaced the hand-rolled binary RAWIV reader/writer in
  `src/vol/RAWIV.cpp` with calls to `cvc::volume`,
  `cvc::volume_file_info`, `cvc::bounding_box`, and `cvc::dimension`
  from libcvc (v3.1.0). Endian handling, header packing, and per-voxel
  type marshalling now flow through libcvc's volume I/O stack.
- Public API of `inc/vol/RAWIV.h` (writeGrid, readRAWIVHeader,
  readShapeCompGrid, readElecGrid) is preserved, so existing F2Dock
  callers (e.g. `Docking.cpp`) compile unchanged.
- The legacy `RAWIVHeader` struct, `bigEndian()` helper, and SWAP
  macros were removed from the header — they are no longer needed.
- `src/vol/CMakeLists.txt` now links `cvc::cvc` PUBLIC on the `vol`
  target so consumers transitively get libcvc's include paths.

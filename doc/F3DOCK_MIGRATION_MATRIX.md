# F3Dock to F2Dock Migration Matrix

This matrix captures a second-pass triage of F3Dock capabilities for incremental porting into modern F2Dock.

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

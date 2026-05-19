# F3Dock to F2Dock Migration Matrix

This matrix captures incremental porting of F3Dock capabilities into modern F2Dock. It is the source of truth for what has landed, what is in progress, and what is deferred.

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

The active mode is printed at startup (`Docking mode: F2Dock (rigid)` / `Docking mode: F3Dock (flexible)`). Defining and parsing the value lives in [inc/f3dock/DockMode.h](../inc/f3dock/DockMode.h); unit tests in [tests/unit/test_dock_mode.cpp](../tests/unit/test_dock_mode.cpp).

All F3Dock-derived stages (`f3dock-loop-closure`, `f3dock-flex`, `f3dock-icp`, `f3dock-domain`) are gated to run **only** when `dockMode == kF3Dock`. Existing F2Dock workflows behave identically with no input changes.

## Capability matrix — pending work

Completed capabilities have moved to the [Completed work](#completed-work) section below. The matrix here lists only items that still need engineering effort.

| Capability | Scientific Value | Effort | Risk | Recommended Action |
|---|---|---:|---:|---|
| Phase 3 follow-ups: real `.f2d` integration test, B-side gridding hoist | Medium | Low-Medium | Low | **Phase 3 close-out (landed)**: real-complex flex integration test landed via PR #22 (`tests/integration/inputs/1A2K_flex_smoke.inp` + baseline, 2-state hinge over `1A2K`); flex-state loop relocated from the CLI into `dock()` / `scoreUntransformed()` (`runFlexStates` in `Docking.cpp`) so library callers also get multi-state behaviour and there is a single seam for a future partition of `dockingMain` into state-invariant setup + per-state body + teardown (tracked as a `TODO(perf)` on `runFlexStates`). |
| Domain graph modeling **runtime integration** | Medium-High for flexible docking workflows | High | High | **Phase 4**: wire `f3dock::domain::DomainGraph` into multi-domain pose composition; needs scoring-loop refactor |
| Blurmaps / GOA molecular surface stack | Medium | High | High | Defer until baseline flexible pipeline lands |
| PRGN / extra geometry stack | Medium | High | High | Defer; extract only needed pieces once use-cases are clear |
| NFFT-based acceleration paths | Medium-High (for some workloads) | High | High | Keep optional; dependency detection already in place, do not hard-require |
| Full F3Dock executable behavior parity | Very High long-term | Very High | Very High | Tracked as the sum of remaining Phase 3 close-out + Phase 4 + deferred items above |

## Porting Principles

1. Keep legacy-derived code optional and isolated behind build options.
2. Add deterministic unit tests before integrating solver outputs into docking flows.
3. Maintain current F2Dock integration regression suite as a release gate.
4. Introduce new scientific dependencies as optional, never mandatory.
5. Use libcvc (currently v3.2.1) as the single source of truth for volume I/O.

## Completed work

### Phase 1 — Tripeptide loop closure (PR #7, merged)

- Added experimental `f3dock-loop-closure` static library.
- Ported tripeptide loop-closure core solver (Coutsias polynomial) into isolated module.
- Stable wrapper API in [inc/f3dock/loop_closure/LoopClosureSolver.h](../inc/f3dock/loop_closure/LoopClosureSolver.h).
- Unit tests for feasible and infeasible closure constraints.
- Optional scientific stack discovery (`Eigen3`, `LAPACK`, `NFFT`).
- DockMode runtime flag (`kF2Dock` / `kF3Dock`) with CLI + input-file parsing.

### Phase 2 — ICP + flex/domain scaffolds (PR #7, merged)

- Added experimental `f3dock-icp` static library (point-to-point ICP).
- Imported `libicp` core components (`matrix`, `kdtree`, point-to-point ICP loop).
- Wrapper API in [inc/f3dock/icp/IcpAligner.h](../inc/f3dock/icp/IcpAligner.h).
- Deterministic unit tests for rigid alignment recovery and input validation.
- Added experimental `f3dock-flex` scaffold for hinge/shear kinematics.
- Deterministic unit tests for hinge rotation and planar shear behavior.
- Added experimental `f3dock-domain` scaffold for domain-graph transform composition.
- Deterministic unit tests for world-transform chaining and cycle rejection.

### Phase 2 completion — libicp parity (PR #7, merged)

- Imported point-to-plane ICP (`libicp::IcpPointToPlane`) including PCA-based normal estimation from the model point cloud.
- Extended `IcpAligner` with `alignPointToPlane(...)` mirroring the point-to-point API, gated by a minimum neighbor count (>= 3).
- Unit tests for surface-fit convergence and input validation (too-few points, too-small neighborhoods).

### Phase 3.5 — libcvc volume I/O adoption (PR #7 + #10, merged)

- Replaced the hand-rolled binary RAWIV reader/writer in [src/vol/RAWIV.cpp](../src/vol/RAWIV.cpp) with calls to `cvc::volume`, `cvc::volume_file_info`, `cvc::bounding_box`, and `cvc::dimension` from libcvc. Endian handling, header packing, and per-voxel type marshalling now flow through libcvc's volume I/O stack.
- Public API of [inc/vol/RAWIV.h](../inc/vol/RAWIV.h) (`writeGrid`, `readRAWIVHeader`, `readShapeCompGrid`, `readElecGrid`) is preserved, so existing F2Dock callers (e.g. `Docking.cpp`) compile unchanged.
- Legacy `RAWIVHeader` struct, `bigEndian()` helper, and SWAP macros removed.
- `src/vol/CMakeLists.txt` links `cvc::cvc` PUBLIC on the `vol` target.
- libcvc bumped to **v3.2.0** (PR #10), picking up SDF / mesher coverage and prebuilt-archive consumption.

### Build / CI hardening (PRs #9, #11, #12, #13, #14, merged)

- Prebuilt libcvc release archives preferred over source build (PR #9).
- Duplicate "Compute artifact metadata" release step removed (PR #11).
- vcpkg cached on Windows; install step retries to absorb gmp/automake mirror flakes (PR #12).
- Full-tree clang-format + CI enforcement (PR #13).
- MSVC C4305 / C4244 narrowing warnings silenced in F2Dock sources (PR #14, + follow-up branch `fix/windows-warnings-2`).

### Phase 3 — hinge/shear runtime integration (PRs #17–#20, merged)

Goal: when `dockMode == kF3Dock`, generate flex-perturbed receptor poses via `f3dock::flex::FlexKinematics` and feed them into the FFT scoring loop. Landed across four stacked PRs:

- **PR #17** — `FlexSamplingConfig` block added to `PARAMS_IN` (hinge axes, angle grid, shear specs). Input-file and CLI parsing wired in; pure-data, no runtime use yet.
- **PR #18** — Startup flex-state enumeration: when `kF3Dock` is active, the `FlexSamplingConfig` is expanded into a deterministic list of `FlexState`s at `dockingMain` entry, gated behind the dock-mode check so F2Dock callers are untouched.
- **PR #19** — Apply flex state 0 (identity) to the receptor in `dockingMain`. Each enumerated state produces a per-state receptor atom array (`f3dock::flex::applyState`); the outer scoring infrastructure consumes the deformed coordinates instead of the raw receptor.
- **PR #20** — Filter init flex-aware via RAII pr-pointer swap. The per-pose forbidden / clash / hydrophobicity filters initialise against the per-state deformed receptor (not the raw rigid receptor) via a scoped pointer swap so the rest of `dockingMain` reads from the active state without invasive parameter threading.

After these four PRs, F3Dock mode runs the rigid scoring pipeline once per flex state with full filter coverage. The outer multi-state loop and per-pose state attribution land in PR #21.

### Phase 3 close-out — outer multi-state loop + per-pose state attribution (PR #21, in review)

- Outer loop in `dockingMain` iterating over every enumerated `FlexState`.
- Per-pose result records tagged with `flex_state_id` so post-filter / re-rank stages can group by state.
- Integration test currently asserts that `kF2Dock` and `kF3Dock` with a trivial single-state config produce identical top hits.

Pending follow-ups tracked in the pending-work matrix above:

- Real `.f2d` complex integration test using `tests/1A2K_L_U.{pqr,f2d}` with a non-identity flex state, asserting the flex top-pose strictly differs from the rigid top-pose.
- Performance optimisation to factor B-side (ligand-side) gridding out of the per-state outer loop: `build_fks` for the ligand, ligand-side `transformAndNormalize`, ligand `gridding` + `griddingElec`, FFT-Kernel construction for B, `getCenterFrequencies` on B, and ligand-side filter object construction all only depend on the rigid ligand inputs (the flex perturbation is applied to the receptor / A-side only). Hoisting them out of `dockingMain`'s per-call work avoids O(N_states - 1) redundant ligand passes.

## Roadmap — remaining phases

### Phase 3 close-out (landed)

Work merged via PR #22 (`feature/f3dock-phase3-followups`) and PR #23 (`feature/f3dock-phase3-flex-loop-in-docking`):

1. **Real-complex integration test** *(landed in PR #22)* — `tests/integration/inputs/1A2K_flex_smoke.inp` exercises `kF3Dock` over the real `1A2K` complex with a z-axis hinge (two non-zero steps -> two flex states). The committed `tests/integration/expected/1A2K_flex_smoke_out.txt` baseline pins the per-state pose tables; the runner uses the existing `compare_f2dock_outputs.py` with `ABS_TOL=1e-3 / REL_TOL=1e-5` over the first 50 rows. Pose records correctly carry the originating state id in the `conf` column (visible in the baseline). Fixes a latent parser bug while at it: `setParamFromFile` was `strtok`'ing on space, so multi-token flex parameters (`flexHingeAxis ax ay az dx dy dz minA stepA n`) only forwarded the first number; the call site now recovers the full value string from the preserved `line` copy before handing it to `FlexSampler::parse_param`.
2. **Flex-state loop relocation into the library entry points** *(landed in PR #23)* — the outer multi-state loop now lives inside `dock()` / `scoreUntransformed()` (`runFlexStates` in `Docking.cpp`) instead of `F2Dock.cpp`'s `main`. The CLI hands off `pr` and the library iterates `pr->flexStates` itself. Behaviour is unchanged for both F2Dock (rigid) and F3Dock (multi-state) modes — verified against the `1A2K`, `1ACB`, and `1A2K_flex` integration baselines. This creates a single seam for the perf hoist below: instead of the CLI invoking `dockingMain` per state, the loop can be split into state-invariant setup + per-state body + teardown without touching any caller.

#### Deferred perf follow-up (open `TODO` on `runFlexStates`)

**B-side gridding hoist** — partition `dockingMain` so that ligand-side setup that is invariant across flex states (`build_fks` for B, FFTW plan creation, B-side gridding + `collectNonzeroGridCells`, kernel construction) runs once per `dock()` call instead of once per state. The rotation loop's per-rotation B work (`transformAndNormalize`, `getCenterFrequencies`, etc.) is not hoistable because it depends on the rotation, not the flex state. Expected saving on N=2 states is modest (a few percent of one state); the larger value of this work is the cleanup of the 1500-line `dockingMain` into three named phases. Tracked as a `TODO(perf, phase3 close-out #2)` comment on `runFlexStates`.

### Phase 4 — domain-graph runtime integration

Goal: replace the implicit single-rigid-body receptor with a multi-domain composition driven by `f3dock::domain::DomainGraph`.

Concrete tasks:
1. ~~Define an input-file syntax for declaring receptor domains and inter-domain joints.~~ **Landed.** `ReceptorDomainConfig`, `DomainSpec`, and `JointSpec` live in [inc/f3dock/domain/DomainSpec.h](../inc/f3dock/domain/DomainSpec.h); `DomainSpecParser::parse_param` accepts three new parameter-file keys (`receptorDomain id firstAtom lastAtom`, `receptorJointFixed parentId childId`, `receptorJointHinge parentId childId ax ay az dx dy dz initAngleRad`), validation (duplicate ids, unknown domains, multi-parent, multi-root) plus `build_graph()` populating a `DomainGraph` with rest-pose joint transforms (Rodrigues rotation about an arbitrary axis line for hinges). `PARAMS_IN::receptorDomains` carries the parsed config; it has no runtime effect yet (atom partitioning + scoring-loop refactor are task 2 / task 3). Unit tests in [tests/unit/test_f3dock_domain_spec.cpp](../tests/unit/test_f3dock_domain_spec.cpp). Diverges from the original "likely a sidecar .domain JSON" suggestion because adding a JSON dependency for this small a schema isn't worth it; the inline parameter-file keys match the existing `flexHingeAxis` / `flexShearPlane` pattern.
2. ~~Replace the receptor's flat atom array with a per-domain partition; compose world transforms via the existing `DomainGraph`.~~ **Landed.** [inc/f3dock/domain/DomainPartition.h](../inc/f3dock/domain/DomainPartition.h) + [src/f3dock-domain/DomainPartition.cpp](../src/f3dock-domain/DomainPartition.cpp) turn the receptor's flat atom array into a dense atom-index → domain-id lookup, validated against the actual receptor size (out-of-range, inverted ranges, and overlapping domains are rejected). `DomainPartition::apply()` consumes a `DomainGraph` whose nodes carry the current per-joint world transforms (rest pose for now; task 3 will swap in per-state poses) and rewrites baseline coordinates by composing each atom's owning world transform; unclaimed atoms ride with the implicit root rigid body unchanged, so receptors without any `receptorDomain` directive are bit-for-bit identical to before. `PARAMS_IN` now carries both the partition and a rest-pose `DomainGraph` alongside the existing `receptorDomains` config; `setParamFromFile` builds and validates them after molecule load so malformed configurations fail at parse time rather than mid-scoring. Unit tests in [tests/unit/test_f3dock_domain_partition.cpp](../tests/unit/test_f3dock_domain_partition.cpp) cover empty / valid / overlap / out-of-range / inverted-range / missing-graph-node scenarios and the per-domain hinge-rotation math.
3. ~~Refactor the FFT scoring loop so each candidate pose is a *vector* of domain transforms rather than a single 6-DoF transform.~~ **Landed.** [inc/f3dock/domain/DomainSampler.h](../inc/f3dock/domain/DomainSampler.h) + [src/f3dock-domain/DomainSampler.cpp](../src/f3dock-domain/DomainSampler.cpp) introduce a `DomainSamplingConfig` of per-hinge sweeps (parsed from the new `receptorJointHingeSweep parentId childId minAngleRad stepRad numSteps` parameter-file key) and a `DomainSampler::enumerate()` that walks the cross-product of hinge angles into a `std::vector<DomainState>` (lexicographic order; last hinge varies fastest; an empty config yields a single empty state so rigid receptors are unaffected). A new `build_graph_with_state()` (in [src/f3dock-domain/DomainSpec.cpp](../src/f3dock-domain/DomainSpec.cpp), declared in DomainSampler.h) re-evaluates the rest-pose graph with each `DomainState`'s hinge-angle overrides, rejecting overrides that name a missing or non-hinge joint. `PARAMS_IN` now carries `domainSampling`, `domainStates`, `domainMaxStates` (default 4096), and `activeDomainStateIndex`; `setParamFromFile` enumerates `domainStates` after the existing flex-state enumeration and prints a one-line summary. The scoring loop `runFlexStates` in [src/f2dock/Docking.cpp](../src/f2dock/Docking.cpp) now iterates the cross-product `n_flex_states × n_domain_states`, setting `activeFlexStateIndex` / `activeDomainStateIndex` per iteration; the two `FlexReceptorGuard*` blocks apply the active domain state via `build_graph_with_state` + `DomainPartition::apply()` on top of any flex-state rewrite, so flex states and domain states compose cleanly (and degenerate to the rigid path when both are empty). Unit tests in [tests/unit/test_f3dock_domain_sampler.cpp](../tests/unit/test_f3dock_domain_sampler.cpp) cover the enumerator (empty / disabled / single / cross-product), the parser (round-trip, unknown keys, bad input), `build_graph_with_state` validation, and an end-to-end partition-with-sampled-state check.
4. ~~Wire `IcpAligner` (point-to-plane) into a refinement pass on the top-N domain-composed poses.~~ **Landed.** [inc/f2dock/IcpRefine.h](../inc/f2dock/IcpRefine.h) + [src/f2dock/IcpRefine.cpp](../src/f2dock/IcpRefine.cpp) expose a `f2dock::refine_pose_point_to_plane(...)` free function that builds a receptor `model` cloud from `xkAOrig/ykAOrig/zkAOrig`, applies the candidate pose's current rotation/translation to the ligand `xkBOrig/ykBOrig/zkBOrig` to form the `moving` cloud, calls `f3dock::icp::IcpAligner::alignPointToPlane(...)` (with neighbor count clamped to `[3, model.size()]` and the same `kMinPoints>=5` guard), and composes the returned ICP delta on the left of the existing pose (`R_new = R_d * R; t_new = R_d * t + t_d`) so refinement is idempotent when ICP returns identity. Four new `PARAMS_IN` fields drive it from the input file (`icpRefine` on/off, `icpRefineTopN` default 10, `icpRefineInlierDist` default 5.0, `icpRefineNumNeighbors` default 10); when `icpRefine` is `false` (the default) the docking path is bit-for-bit unchanged. Inside `printIntermediateStats` in [src/f2dock/Docking.cpp](../src/f2dock/Docking.cpp) the refinement runs only for the first `icpRefineTopN` ranked poses, mutates the in-memory `transformation` matrix before the line is written to the `.out` file, and emits a single `ICP refine rank=... rmsd_before=... rmsd_after=...` line to stdout per refined pose (the `.out` schema is unchanged so existing parsers keep working). Unit tests in [tests/unit/test_f2dock_icp_refine.cpp](../tests/unit/test_f2dock_icp_refine.cpp) cover identity-pose stability, small-translation recovery (z-component shrinks toward 0), too-few-points rejection, ICP-delta composition with an existing rotated pose (with `det(R)≈1` check), and null-pointer rejection.
5. Integration test: a hinge-joined two-domain receptor against a rigid ligand, comparing F3Dock-mode top hit against a pre-computed reference.

Risk: high. This is the scoring-loop refactor that the original F3Dock executable was built around; expect this to span multiple PRs.

### Deferred until after Phase 4

- Blurmaps / GOA molecular surface stack.
- PRGN / extra geometry stack.
- NFFT-accelerated scoring paths (dependency detection is already in place; only enable behind an explicit flag after Phase 4 lands).

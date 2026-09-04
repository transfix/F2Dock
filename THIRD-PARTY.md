# Third-party code in F2Dock

F2Dock is distributed under the GNU Lesser General Public License, version 2.1
only — see [`LICENSE`](LICENSE). Most of the tree is that: work by the
Computational Visualization Center at The University of Texas at Austin,
carrying the UT Austin LGPL-2.1 banner.

The rest is vendored from elsewhere. This file inventories every such tree with
its upstream, author and licence, and says whether the default build compiles
it.

This is a record of what the source files say. It is not legal advice, and two
entries below (`libmol`, `ElementInformation.h`) are open questions rather than
settled facts.

## Summary

| Tree | Files | Licence | In the default build? |
|---|---:|---|---|
| [`inc/libmol`, `src/libmol`](#libmol) | 76 | Boston University notice — **not free, not redistributable** | No (`HAVE_LIBMOL=OFF`) |
| [`inc/fft-utils/sparsefft3*`, `src/fft-utils/sparsefft3*`](#sparsefft3) | 4 | GPL-2.0-or-later | No (`F2DOCK_ENABLE_GPL_COMPONENTS=OFF`) |
| [`inc/libicp`, `src/f3dock-icp`](#libicp) | 8 | GPL-2.0-or-later, plus AFL-1.1 | No (`F2DOCK_ENABLE_GPL_COMPONENTS=OFF`) |
| [`src/f3dock-loop-closure/tripep_closure.*`](#tripeptide-loop-closure) | 2 | **No grant stated** | Yes |
| [`src/f3dock-loop-closure/sturm.h`](#sturm-solver) | 1 | Graphics Gems, grant not restated | Yes |
| [`inc/fast-GB/SSEApproxMath.*`](#sse-approximate-math) | 2 | zlib | Yes |
| [`inc/fft-utils/fftw3.h`](#vendored-fftw3h) | 1 | BSD-2-Clause (header only) | Yes |
| [`inc/PG-range/cuckoo.h`, `src/PG-range/cuckoo.cc`](#cuckoo-hashing) | 2 | Permissive, attribution required | Yes |
| [`inc/PG-range/dsexceptions.h`](#dsexceptionsh) | 1 | **Provenance unknown** | Yes |

External libraries F2Dock links but does not vendor — FFTW3, pthreads,
libcvc, and the optional Eigen3 / LAPACK / NFFT3 — are not listed here. FFTW3
is itself GPL-2.0-or-later and is linked dynamically.

## The GPL build switch

Two of the vendored trees are GPL-2.0-or-later. Linking either into F2Dock
makes the resulting binary a combined work covered by the GPL, which is not
what this project distributes, so both are behind one option that is **OFF by
default**:

```
cmake -S . -B build -DF2DOCK_ENABLE_GPL_COMPONENTS=ON
```

With the option OFF the GPL sources are not compiled and not linked, and:

- `useSparseFFT` defaults to `false` and is refused if a parameter file asks
  for it. Transforms go through dense FFTW instead — slower, same scores.
- `icpRefine` is refused if a parameter file asks for it, and
  `f2dock::refine_pose_point_to_plane` is a no-op returning `refined == false`.

With the option ON, CMake emits a warning and the binary reports
`build licence: GPL-2.0-or-later` at startup. Do not ship that binary under
F2Dock's own licence.

The LGPL core does not `#include` the GPL headers in either configuration.
`inc/fft-utils/sparsefft3-api.h` is the seam: it forwards to the real header
when the module is built and declares the same entry points against an
incomplete type otherwise.

---

## libmol

- **Path:** `inc/libmol/` (38 headers), `src/libmol/` (38 sources + CMake)
- **Upstream:** Structural Bioinformatics Laboratory, Boston University
- **Contact of record:** Dima Kozakov, Department of Biomedical Engineering,
  `midas@bu.edu`
- **Licence:** none of the usual ones. 72 of the 76 files carry the notice at
  `inc/libmol/mask.h:3`–`16`; the remaining four (`newgetline.{h,cpp}`,
  `test.cpp`, `testHbond.cpp`) have no header, but they are part of the same
  vendored tree.

The notice grants copy-and-modify rights **only** to "US Universities and US
Government supported Research Institutions", **only** for "educational and
research purposes", and states that the software "may not be distributed to
any other institution without express permission from the Structural
Bioinformatics Laboratory and Boston University".

That is a field-of-use and class-of-licensee restriction. It is incompatible
with the LGPL, incompatible with the GPL, and incompatible with redistribution
generally.

`HAVE_LIBMOL` is OFF by default (`src/CMakeLists.txt`), so the tree is not
compiled — but it is present in the published source tree, and source
redistribution is exactly what the notice restricts. Neither `CVC-Lab/F2Dock`
nor `CVC-Lab/F3Dock` ships `libmol`, while `CVC-Lab/F3Dock`'s
`src/CMakeLists.txt` still calls `ADD_SUBDIRECTORY (libmol)` for a directory
that is not in its checkout — consistent with the tree having been stripped
before those repositories were published.

`src/hbondFilter/` is UT Austin's own LGPL code, but it only builds when
`HAVE_LIBMOL` is ON because it consumes libmol.

**Status: open.** Removing the tree, or obtaining permission for it, is a
decision for UT — not something to settle in a commit. Until it is settled,
this repository should be treated as not redistributable to third parties.

## sparsefft3

- **Path:** `inc/fft-utils/sparsefft3.h`, `inc/fft-utils/sparsefft3-timers.h`,
  `src/fft-utils/sparsefft3.cpp`, `src/fft-utils/sparsefft3-plan.cpp`
- **Upstream:** Massachusetts Institute of Technology, 2000 (the sparse 3-D
  FFT distributed alongside FFTW); modified at CVC Lab, UT Austin, for FFTW3
- **Licence:** GPL-2.0-or-later (`inc/fft-utils/sparsefft3.h:5`–`9`)
- **Used by:** the `useSparseFFT` parameter — the sparse forward/backward
  transform pair over the frequency maps

Not built by default. See [the GPL build switch](#the-gpl-build-switch).

Replacing this with a permissively-licensed sparse transform is the cleanest
long-term fix, since the dense fallback is slower.

## libicp

- **Path:** `inc/libicp/` (5 headers), `src/f3dock-icp/` (6 sources)
- **Upstream:** Institute of Measurement and Control Systems, Karlsruhe
  Institute of Technology
- **Author:** Andreas Geiger
- **Licence:** GPL-2.0-or-later (`inc/libicp/icp.h:8`–`10`) for `icp`,
  `icpPointToPoint`, `icpPointToPlane` and `matrix`
- **Used by:** the `icpRefine` parameter — point-to-plane refinement of the
  top-N final poses

`inc/libicp/kdtree.h` and `src/f3dock-icp/kdtree.cpp` are different: they are
(c) Matthew B. Kennel, Institute for Nonlinear Science, UCSD (2004), under the
**Academic Free License version 1.1**, and the header points at "file LICENSE
with additional provisions in that same file" — **a file that is not in this
tree**. The additional provisions are therefore unknown. AFL-1.1 is also not
GPL-compatible, so `kdtree` sitting inside a GPL library is its own problem.

`src/f3dock-icp/IcpAligner.cpp` and `inc/f3dock/icp/IcpAligner.h` are UT
Austin's own LGPL wrapper, but they have nothing to wrap without libicp, so
the whole target is gated together.

Not built by default. See [the GPL build switch](#the-gpl-build-switch).

## Tripeptide loop closure

- **Path:** `src/f3dock-loop-closure/tripep_closure.h`,
  `src/f3dock-loop-closure/tripep_closure.cpp`
- **Authors:** Chaok Seok, Evangelos Coutsias, Matthew Jacobson, Ken Dill
  (UCSF and University of New Mexico), 2003; written by Chaok Seok 2003, ported
  to C++ by Nathan Clement 2017
- **Licence:** **none stated.** The files carry a bare copyright line and no
  grant of any kind.

A bare copyright notice with no licence grants no rights. This is built by
default, via `F2DOCK_ENABLE_F3DOCK_EXPERIMENTAL`.

**Status: open.** The algorithm is published (Coutsias et al., *J. Comput.
Chem.* 25:510, 2004) and the original Fortran/C is widely redistributed, which
suggests the omission is an oversight rather than a restriction — but that is a
guess, and the fix is to ask the authors for an explicit grant.

## Sturm solver

- **Path:** `src/f3dock-loop-closure/sturm.h`, and the `sturm` routines
  compiled with the loop-closure module
- **Upstream:** "Using Sturm Sequences to Bracket Real Roots of Polynomial
  Equations", D. G. Hook and P. R. McAree, *Graphics Gems*, Academic Press,
  1990, via `acm.org/pubs/tog/GraphicsGems`; modified by Chaok Seok 2003 and
  Nathan Clement 2017
- **Licence:** not restated in this tree. The Graphics Gems series
  distributed its code under a blanket permission notice from the publisher,
  but that notice does not appear here, and these copies have been modified
  twice since.

Built by default. Left without a header deliberately: stamping F2Dock's own
banner on it would misattribute it.

## SSE approximate math

- **Path:** `inc/fast-GB/SSEApproxMath.h`, `src/fast-GB/SSEApproxMath.cpp`
- **Author:** Julien Pommier, 2007 (SIMD `exp`/`log`, after the Intel
  Approximate Math library and cephes)
- **Licence:** zlib

Permissive, attribution preserved in the file, LGPL-compatible. Built by
default. No action needed.

## Vendored fftw3.h

- **Path:** `inc/fft-utils/fftw3.h`
- **Authors:** Matteo Frigo and MIT, 2003/2006
- **Licence:** BSD-2-Clause. The file says so explicitly: the BSD terms apply
  "*only* to this header file, and *not* to the other files distributed with
  FFTW or derived therefrom". The FFTW library itself is GPL-2.0-or-later and
  is linked dynamically.

Permissive and LGPL-compatible. It does shadow the system `fftw3.h`, which is a
maintenance hazard independent of licensing.

### The library behind that header is the largest open licensing question

The *header* is BSD-2. **The library is GPL-2.0-or-later**, and the cvcpkg
recipe agrees (`recipes/fftw3/recipe.yaml`: `license: GPL-2.0-or-later`).
`SetupFFTW()` is called unconditionally at `CMakeLists.txt:32` and `FFTW_LIB`
is linked unconditionally, so **every F2Dock binary is a combined work with
GPL code** — independently of the vendored `sparsefft3` and `libicp` this file
already covers, and not fixed by gating those off. Under the usual reading,
dynamic linking does not change that.

That makes an FFT-provider swap the load-bearing item for LGPL status, not an
optimisation. The swap itself is small: F2Dock's whole FFTW surface is already
funnelled through the `FFTW_*` macros in `inc/fft-utils/fftwPrecision.h`, and
amounts to `plan_dft_3d` (11 sites), `plan_dft_1d` (5), `plan_dft_r2c_3d` (4),
`plan_dft_c2r_3d` (3), `plan_dft_2d` (2), plus `execute`/`destroy_plan`/
`malloc`/`free`. No threads, no guru interface, no wisdom, no advanced plans.

**Measured, so the trade is a number rather than a guess.** PocketFFT
(BSD-3-Clause, the transform NumPy and SciPy use) against FFTW 3, double
precision, 3-D complex-to-complex, plan-once/execute-many, on this workstation:

| N | FFTW s/exec | PocketFFT s/exec | ratio | max abs diff |
|---|---|---|---|---|
| 64 | 0.00295 | 0.00435 | 1.47x | 5.8e-13 |
| 100 | 0.01507 | 0.02280 | 1.51x | 1.6e-12 |
| **120** (F2Dock's real `numFreq`) | **0.03535** | **0.04640** | **1.31x** | 2.3e-12 |
| 128 | 0.03444 | 0.06339 | 1.84x | 2.2e-12 |
| 200 | 0.17976 | 0.22878 | 1.27x | 4.9e-12 |
| 256 | 0.36164 | 0.59289 | 1.64x | 7.3e-12 |

Numerically the two agree to float64 round-off, so PocketFFT is a valid drop-in
on correctness. On speed, the honest figure is *throughput at equal core
count*, because F2Dock parallelises over rotations with single-threaded
transforms rather than threading inside one. Measuring W concurrent
single-threaded transforms at N=120:

| workers | FFTW (xf/s) | PocketFFT (xf/s) | FFTW faster by |
|---|---|---|---|
| 1 | 26.9 | 16.2 | 1.66x |
| 2 | 53.2 | 30.7 | 1.73x |
| 4 | 102.1 | 67.9 | 1.50x |
| 8 | 161.3 | 114.0 | 1.41x |

So **PocketFFT costs roughly 1.4–1.7x FFT throughput** in F2Dock's actual
architecture. (A single-transform latency comparison flatters PocketFFT —
`nthreads=4` makes one transform 2.6x *faster* than single-threaded FFTW — but
that is spending four cores to do one transform, and F2Dock already spends
those cores on four rotations.)

Options, in order of licensing cleanliness:

1. **PocketFFT (BSD-3)** — clean LGPL, header-only, no new bundle needed, at
   ~1.4–1.7x FFT cost.
2. **Intel MKL or cuFFT** — proprietary but redistributable, which LGPL-2.1
   permits linking; likely matches or beats FFTW. This is the cvcpkg roadmap's
   D9 item.
3. **Buy FFTW's commercial licence from MIT** — no code change at all.

`sparsefft3` sits on top of whichever base transform is chosen: it is a pruned
FFT exploiting input/output sparsity, and its output is *identical* to the
dense path (verified on 1A2K), so it is recoverable as a pure optimisation
against a permissive base rather than something that must be kept.

## Cuckoo hashing

- **Path:** `inc/PG-range/cuckoo.h`, `src/PG-range/cuckoo.cc`
- **Authors:** Rasmus Pagh and Flemming Friche Rodler, BRICS, Department of
  Computer Science, University of Aarhus, 2001
- **Licence:** permissive — use, copy, modify and distribute for any purpose
  without fee, provided due acknowledgement to the authors is given and the
  permission notice appears in all copies.

LGPL-compatible. Built by default. The acknowledgement requirement is met by
this file and by the notice in the sources.

## dsexceptions.h

- **Path:** `inc/PG-range/dsexceptions.h`
- **Licence:** none stated; nine lines declaring five empty exception classes.

Arrived with the original 2012 import. The name and content match a
long-circulating teaching header, so it is probably not original to this
project — hence no F2Dock banner on it. Given that it is five empty class
declarations there is likely nothing copyrightable in it either, but the
provenance is genuinely unknown.

**Status: minor, open.** The cheap fix is to replace it with an equivalent
written here.

---

## Two things this file does not settle

**`inc/f2dock/ElementInformation.h` provenance.** The same file is
GPL-2.0-or-later in both `CVC-Lab/F2Dock` and `CVC-Lab/F3Dock`, and
LGPL-2.1-only here, with the copyright year changed 2003 → 2011. UT Austin
holds both, so UT can legitimately relicense its own work — but nothing in the
tree or in the git history records that decision. It needs a note from whoever
made it, not a commit from whoever noticed.

**`libmol`.** See above. It is the sharpest item in this file.

# F2Dock Overview

## What is F2Dock?

F2Dock is a **rigid-body protein-protein docking program** developed at the University of Texas at Austin's Computational Visualization Center. It predicts binding poses between two proteins by systematically exploring rotational and translational configurations in 3D space, using FFT-accelerated scoring to efficiently evaluate millions of candidate poses.

The distribution also includes **GB-Rerank**, a companion tool that rescores the top predicted docking poses using the more computationally expensive Generalized Born (GB) electrostatic solvation model.

**Authors:** Rezaul Alam Chowdhury, Muhibur Rasheed (advisor: Chandrajit Bajaj)  
**License:** GNU Lesser General Public License v2.1  
**Usage:**
```
F2Dock -score|saveGrid|vdw|effGridFile parameterFile
GB-rerank parameterFile
```

---

## Project State

- **Build system:** CMake 2.6+ (legacy, but functional)
- **Language:** C++ with some C utilities
- **Threading:** POSIX pthreads
- **Known issues:**
  - `CMakeLists.txt` has a hardcoded include path (`/User/andre/F2Dock-refactored/inc/`) that needs fixing
  - `CMakeCache.txt` is stale (from a macOS build at a different path)
  - Build artifacts (Makefiles, CMakeFiles/) are committed to the repo
  - Vim swap files (`.cpp.swp`) present in source tree
- **Overall:** Source structure is complete and well-organized. Needs minor CMake fixes for a clean build on a modern Linux system.

---

## Entry Points / Executables

| Executable       | Source File                          | Purpose                                          |
|------------------|--------------------------------------|--------------------------------------------------|
| **F2Dock**       | `src/f2dock/F2Dock.cpp`              | Main docking engine (FFT-accelerated scoring)    |
| **GB-Rerank**    | `src/GB-rerank/GB-rerank.cpp`        | Rescores docking results using GB solvation      |
| **F2DockServer** | `src/F2DockServer/F2DockServer.cpp`  | XML-RPC server for remote docking job submission |
| **FilterOutput** | `src/PostFiltering/FilterOutput.cpp` | Post-processing filter for docking results       |

---

## Source Structure

### Core Infrastructure
| Module        | Directory          | Purpose                                    |
|---------------|--------------------|--------------------------------------------|
| **math**      | `src/math/`        | Linear algebra: Matrix, Vector, Quaternion, Gaussian, Ray |
| **fft-utils** | `src/fft-utils/`   | FFTW wrappers and custom sparse FFT        |
| **utils**     | `src/utils/`       | String/file handling, memory management    |
| **vol**       | `src/vol/`         | Volume data I/O (RAWIV format)             |
| **misc-ident**| `src/misc-ident/`  | Element/atom property lookup               |

### Scoring & Filtering Modules
| Module          | Directory             | Purpose                                   |
|-----------------|-----------------------|-------------------------------------------|
| **fast-clash**  | `src/fast-clash/`     | Van der Waals clash detection             |
| **fast-GB**     | `src/fast-GB/`        | Generalized Born electrostatic solvation  |
| **fast-hydro**  | `src/fast-hydro/`     | Hydration free energy scoring             |
| **fast-LJ**     | `src/fast-LJ/`       | Lennard-Jones potential scoring           |
| **fast-PQ**     | `src/fast-PQ/`       | Partial charge electrostatics             |
| **fast-resCont**| `src/fast-resCont/`   | Residue contact filtering                 |

### Search & Ranking
| Module        | Directory          | Purpose                                    |
|---------------|--------------------|--------------------------------------------|
| **PG-range**  | `src/PG-range/`    | Spatial indexing / geometric range queries  |

### Optional Modules (compile-time flag `HAVE_LIBMOL`)
| Module          | Directory            | Purpose                                  |
|-----------------|----------------------|------------------------------------------|
| **libmol**      | `src/libmol/`       | Molecular structure library              |
| **hbondFilter** | `src/hbondFilter/`  | Hydrogen bond detection and filtering    |

### Networking & Server
| Module          | Directory              | Purpose                                |
|-----------------|------------------------|----------------------------------------|
| **XmlRPC**      | `src/XmlRPC/`         | XML-RPC protocol library               |
| **F2DockServer**| `src/F2DockServer/`   | Remote docking job server              |

### Post-Processing
| Module            | Directory                | Purpose                             |
|-------------------|--------------------------|-------------------------------------|
| **PostFiltering** | `src/PostFiltering/`     | Result filtering and formatting     |

---

## API / Key Classes

Headers live under `inc/` with one subdirectory per module:

- **`inc/f2dock/`** — `Docking.h`, `TopValues.h`, `ValuePosition3D.h`, `ElementInformation.h`
- **`inc/math/`** — `Matrix.h`, `Vector.h`, `Quaternion.h`, `Gaussian.h`, `SmoothingFunction.h`
- **`inc/fft-utils/`** — `fftw3.h`, `fftwPrecision.h`, `rank-fftw.h`, `sparsefft3.h`, `fastfft.h`
- **`inc/fast-GB/`** — `fastGpol.h`, `fastBornRadius.h`, `fastDispE.h`, `SSEApproxMath.h`
- **`inc/fast-clash/`** — `clashFilter.h`
- **`inc/vol/`** — `RAWIV.h`
- **`inc/utils/`** — General utilities
- **`inc/XmlRPC/`** — XML-RPC client/server headers

Key classes: `Docking` (main orchestrator), `Matrix`/`Vector`/`Quaternion` (geometry), `fastGpol`/`fastBornRadius` (GB solvation), `clashFilter` (clash detection).

---

## Dependencies

### Required
- **FFTW3** — Fast Fourier Transform library (double-precision preferred; single-precision fallback)
  - Optional threading support via `fftw3_threads`
- **pthreads** — POSIX threading
- **libc, libm** — Standard C runtime and math

### Optional
- **libmol** — Enables hydrogen bond filtering (`-DHAVE_LIBMOL=ON`)

### Build System
- **CMake** ≥ 2.6
- **C++ compiler** with C++03 support (GCC, Clang)

---

## Tests

Four protein complex test cases in `tests/`:

| PDB ID | Files Provided |
|--------|----------------|
| **1A2K** | `.f2d`, `.pqr`, `.pqr.pdb`, `.pqr.rawn`, `.quad`, `.inp`, `_out.txt`, `_rmsd_backbone_unbound_10.0.txt` |
| **1ACB** | Same set |
| **1AHW** | Same set (no `_out.txt`) |
| **1JZD** | Same set (no `_out.txt`) |

Each `.inp` file is a docking configuration specifying receptor/ligand files, rotation sampling, scoring parameters, and thread count. The `_out.txt` files contain expected output for validation. The `_rmsd_backbone_unbound_10.0.txt` files contain RMSD data for evaluating prediction quality.

Example (from `1A2K.inp`): 100,000 rotations, 20,000 output poses, GB reranking of top 2,000, 4 threads.

**Note:** These are integration-level tests (run the full docking pipeline on known complexes), not unit tests. There is no automated test harness — tests are run manually by comparing output against `_out.txt`.

---

## Build Instructions

```bash
cd /home/joe/src/cvc/F2Dock
mkdir -p build && cd build
cmake ..
make -j$(nproc)
```

Executables will be placed in `build/bin/`, libraries in `build/lib/`.

---

## Publications

References from the [CVC Lab F2Dock project page](https://cvc-lab.github.io/software/f2dock/) and the [protein docking project page](http://www.cs.utexas.edu/~bajaj/cvc/projects/angstrom/sd/dock.shtml):

### Primary F2Dock Papers

1. C. Bajaj, R. Chowdhury, and V. Siddavanahalli.
   **F2Dock: Fast Fourier Protein-Protein Docking.**
   *IEEE/ACM Transactions on Computational Biology and Bioinformatics*, 8(1):45–58, 2011.
   [DOI](http://doi.ieeecomputersociety.org/10.1109/TCBB.2009.57)

2. R. Chowdhury, D. Keidel, M. Moussalem, A. Olson, M. Sanner, and C. Bajaj.
   **F2Dock 2.0: Improved Fast Fourier Protein-Protein Docking.**
   Under Preparation, 2010.

3. D. Keidel, R. Chowdhury, M. Moussalem, A. Olson, M. Sanner, and C. Bajaj.
   **F2Dock 2.0 Scoring: An Evaluation.**
   Under Preparation, 2010.

### Flexible Docking

4. C. Bajaj, R. Chowdhury, and V. Siddavanahalli.
   **F3Dock: A Fast, Flexible and Fourier Based Approach to Protein-Protein Docking.**
   The University of Texas at Austin, ICES Report 08-01, January 2008.

### Reranking / Solvation Energy

5. M. Moussalem, R. Chowdhury, D. Keidel, A. Olson, M. Sanner, and C. Bajaj.
   **Comparative Study of Three Reranking Methods for Fast Fourier Transform-Based Protein-Protein Docking Programs.**
   Under Preparation, 2010.

### Supporting Methods

6. C. Bajaj, R. Chowdhury, and M. Rasheed.
   **A Dynamic Data Structure for Flexible Molecular Maintenance and Informatics.**
   *Proceedings of the ACM Symposium on Solid and Physical Modeling (SPM 2009)*, San Francisco, California, pp. 259–270, 2009.
   [ACM DL](http://portal.acm.org/citation.cfm?id=1248377.1248392)

7. C. Bajaj, S.-C. Chen, and A. Rand.
   **An Efficient Higher-Order Fast Multipole Boundary Element Solution for Poisson-Boltzmann Based Molecular Electrostatics.**
   *SIAM Journal on Scientific Computing*, 2011.

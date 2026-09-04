# F2Dock

[![CI/CD](https://github.com/transfix/F2Dock/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/transfix/F2Dock/actions/workflows/ci.yml)

Fast Fourier protein–protein docking with multi-term scoring.

F2Dock predicts how two proteins bind by exhaustively sampling rigid-body orientations on an FFT grid and ranking poses with a combined electrostatic, van der Waals, shape complementarity, and solvation score.  A companion tool, **GB-Rerank**, rescores the top poses using a Generalized Born solvation model for improved accuracy.

Developed at the [Computational Visualization Center](https://cvc-lab.github.io), Oden Institute for Computational Engineering and Sciences, The University of Texas at Austin.

## Features

- FFT-accelerated exhaustive rotational/translational search
- Multi-term scoring: electrostatics, van der Waals, desolvation, shape complementarity
- Generalized Born reranking (GB-Rerank)
- Hydrogen-bond and residue-contact post-filters
- Dynamic Packing Grid (PG) for fast spatial range queries
- Cross-platform: Linux, macOS, Windows

## Quick start

### Prerequisites

| Dependency | Notes |
|---|---|
| C++20 compiler | GCC 13+, Clang 16+, or MSVC 19.4+ |
| Python 3.9+ | Only to run `cvcpkg` itself |
| [cvcpkg](https://cvcpkg.org) | `pip install cvcpkg` |

Everything else — libcvc, FFTW3, Boost, CGAL, HDF5, ImageMagick, pthreads4w
on Windows, and even CMake and Ninja — is installed from the cvcpkg catalog.
F2Dock does not use apt, Homebrew or vcpkg, and does not build libcvc from
source.

### Build

```bash
pip install cvcpkg
cvcpkg install-deps cvcpkg/recipes/f2dock --prefix deps \
    --config release --link shared --include-host-tools
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$PWD/deps"
cmake --build build --parallel
```

Executables are placed in `build/bin/`. The dependency set is whatever
[`cvcpkg/recipes/f2dock/recipe.yaml`](cvcpkg/recipes/f2dock/recipe.yaml)
declares, resolved transitively — the same recipe CI uses, and the same one
that builds and publishes the `f2dock` package.

### Packaging

F2Dock ships as a cvcpkg package under the `cvc` organization:

```bash
cvcpkg install f2dock --prefix /some/prefix --config release --link shared
```

#### Optional GPL components

Two vendored modules are GPL-2.0-or-later while F2Dock itself is LGPL-2.1-only,
so they are **not built by default**:

| Module | Feature it backs | Upstream |
|---|---|---|
| `sparsefft3` | `useSparseFFT` — sparse 3-D FFT | MIT 2000, modified at CVC Lab |
| `libicp` | `icpRefine` — point-to-plane pose refinement | Andreas Geiger, KIT |

The default build uses the **dense FFTW** transform instead of the sparse one.
It is slower, and it produces the same scores. `icpRefine` is unavailable.
Both parameters are rejected with an explanatory message rather than silently
ignored, and the binary prints which modules it contains at startup:

```
GPL components: sparseFFT=no libicp=no (build licence: LGPL-2.1-only)
FFT path: dense (FFTW)
```

To build them:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DF2DOCK_ENABLE_GPL_COMPONENTS=ON
```

The resulting binary is a combined work covered by the GNU GPL v2 or later,
**not** by F2Dock's LGPL-2.1. See [THIRD-PARTY.md](THIRD-PARTY.md).

### Run

```
F2Dock parameterFile
F2Dock -score parameterFile
GB-Rerank parameterFile
```

## Usage example

F2Dock reads a plain-text parameter file that lists the input molecules and scoring options.  Four example parameter files and their associated test data are included in `tests/`.

Run docking on the 1ACB enzyme–inhibitor complex from the repository root:

```bash
build/bin/F2Dock tests/1ACB.inp
```

The parameter file `tests/1ACB.inp` looks like:

```
staticMolecule   tests/1ACB_R_U_1.7.f2d
movingMolecule   tests/1ACB_L_U.f2d
staticMoleculePQR  tests/1ACB_R_U.pqr
movingMoleculePQR  tests/1ACB_L_U.pqr
staticMoleculeQUAD tests/1ACB_R_U.quad
movingMoleculeQUAD tests/1ACB_L_U.quad
rmsdAtoms        tests/1ACB_rmsd_backbone_unbound_10.0.txt
outFile          tests/1ACB_out.txt

rotFile          deg15.mtx
numRot           100000
effGridFile      fftw-rank.txt

complexType      E
numSolutions     20000
rerank           true
numRerank        2000
breakDownScores  true
numThreads       4
```

### Key parameters

| Parameter | Description |
|---|---|
| `staticMolecule` / `movingMolecule` | F2D structure files for receptor and ligand |
| `staticMoleculePQR` / `movingMoleculePQR` | PQR charge/radius files |
| `staticMoleculeQUAD` / `movingMoleculeQUAD` | Surface quadrature point files |
| `rotFile` | Euler angle rotation set (e.g. `deg15.mtx`) |
| `numRot` | Number of rotations to sample |
| `effGridFile` | Machine-specific FFT grid sizes for best performance |
| `complexType` | `A` (antibody–antigen), `E` (enzyme–inhibitor), `G` (other), `U` (auto-detect) |
| `numSolutions` | Max docking poses to output |
| `rerank` / `numRerank` | Enable GB-Rerank rescoring and number of top poses to rescore |
| `outFile` | Output file for ranked docking poses |
| `numThreads` | Number of CPU threads |
| `useSparseFFT` | Sparse FFT instead of dense FFTW. Requires `-DF2DOCK_ENABLE_GPL_COMPONENTS=ON`; see [Optional GPL components](#optional-gpl-components) |
| `icpRefine` | Point-to-plane ICP refinement of the top poses. Requires `-DF2DOCK_ENABLE_GPL_COMPONENTS=ON` |

All file paths are relative to the current working directory.  See [doc/params/](doc/params/) for the full parameter reference.

## Windows

```bash
vcpkg install fftw3:x64-windows pthreads:x64-windows
cmake -B build -DCMAKE_TOOLCHAIN_FILE=%VCPKG_ROOT%/scripts/buildsystems/vcpkg.cmake
cmake --build build --config Release
```

## Testing

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
cd build && ctest --output-on-failure
```

89 unit tests cover the spatial index, scoring components, file parsers, and math utilities.

## Formatting

F2Dock CI enforces clang-format on changed lines.

Use this before pushing:

```bash
git fetch origin
git clang-format <base>
```

For example, if your PR targets `master`:

```bash
git clang-format origin/master
```

## Install

```bash
cmake --install build --prefix /usr/local
```

Installs executables to `bin/` and data files (rotation matrices) to `share/F2Dock/data/`.

## Releases

Pre-built binaries for Linux, macOS, and Windows are published on the [Releases](https://github.com/transfix/F2Dock/releases) page whenever a version tag is pushed.

## Project layout

```
src/
  f2dock/         Main docking driver
  fast-GB/        Fast Generalized Born solvation
  fast-LJ/        Lennard-Jones van der Waals
  fast-hydro/     Hydrophobicity scoring
  fast-PQ/        Skin surface / shape complementarity
  fast-clash/     Steric clash filter
  fast-resCont/   Residue contact filter
  PG-range/       Dynamic Packing Grid spatial index
  fft-utils/      FFT grid construction and correlation
  GB-rerank/      GB-Rerank standalone rescoring tool
  hbondFilter/    Hydrogen-bond post-filter
  math/           Pairing heap, spherical harmonics
  libmol/         Molecule I/O (PDB/PQR/F2D) -- not built; see THIRD-PARTY.md
  utils/          Parameter parsing, timers
inc/              Headers (mirrors src/ structure)
tests/unit/       Google Test suite
doc/              Parameter and dependency documentation
```

## References

- Bajaj, C., Chowdhury, R., and Rasheed, M. **A Dynamic Data Structure for Flexible Molecular Maintenance and Informatics.** *Bioinformatics*, 27(1):55–62, 2011.
- Chowdhury, R., Rasheed, M., Keidel, D., Moussalem, M., Olson, A., Sanner, M., and Bajaj, C. **Protein-Protein Docking with F2Dock 2.0 and GB-Rerank.** *PLoS ONE*, 8(3):e51397, 2013.

## License

[GNU Lesser General Public License v2.1](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
only. The full text is in [`LICENSE`](LICENSE), and it is what CPack ships in
generated packages.

Vendored third-party code, its licences, and which of it the default build
compiles are inventoried in [THIRD-PARTY.md](THIRD-PARTY.md). Two entries there
are unresolved and are flagged as such.

## Contact

[Computational Visualization Center](https://cvc-lab.github.io)
Oden Institute for Computational Engineering and Sciences
The University of Texas at Austin

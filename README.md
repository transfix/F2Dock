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
| CMake 3.16+ | |
| FFTW 3 | `libfftw3-dev` (apt), `fftw` (brew), or via vcpkg |
| pthreads | Built-in on Linux/macOS; PThreads4W via vcpkg on Windows |

### Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

Executables are placed in `build/bin/`.

### Run

```
F2Dock --help
F2Dock parameterFile
F2Dock -score parameterFile
GB-Rerank parameterFile
```

See the [doc/params/](doc/params/) directory for parameter file documentation.

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
  libmol/         Molecule I/O (PDB/PQR/F2D)
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

## Contact

[Computational Visualization Center](https://cvc-lab.github.io)
Oden Institute for Computational Engineering and Sciences
The University of Texas at Austin

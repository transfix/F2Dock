#!/usr/bin/env bash
# cvcpkg/recipes/f2dock/build.sh — build F2Dock against an installed libcvc SDK.
#
# libcvc is NOT rebuilt here: it comes from the deps prefix as a published
# cvcpkg bundle, which is the whole point. Building it from source is what
# dragged CGAL and the platform toolchains into F2Dock's build and broke
# macOS (CGAL 6.2) and Windows (v180 CL.exe) -- see the top-level
# CMakeLists.txt comment on find_package(cvc).
set -euo pipefail

: "${CVC_SOURCE_DIR:?CVC_SOURCE_DIR must be set}"
: "${CVC_BUILD_DIR:?CVC_BUILD_DIR must be set}"
: "${CVC_INSTALL_DIR:?CVC_INSTALL_DIR must be set}"

CVC_BUILD_TYPE="${CVC_BUILD_TYPE:-Release}"
CVC_JOBS="${CVC_JOBS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)}"

case "$(echo "$CVC_BUILD_TYPE" | tr '[:upper:]' '[:lower:]')" in
  debug) CMAKE_BUILD_TYPE=Debug ;;
  *) CMAKE_BUILD_TYPE=Release ;;
esac

CMAKE_ARGS=(
  -G Ninja
  -S "$CVC_SOURCE_DIR"
  -B "$CVC_BUILD_DIR"
  -DCMAKE_INSTALL_PREFIX="$CVC_INSTALL_DIR"
  -DCMAKE_BUILD_TYPE="$CMAKE_BUILD_TYPE"
  # Tests pull googletest over the network; a recipe build is not the place
  # for that. CI runs the suite separately from a normal source build.
  -DF2DOCK_BUILD_TESTS=OFF
  # The GPL modules (sparsefft3, libicp) stay OFF so the published bundle is
  # unambiguously LGPL-2.1. See THIRD-PARTY.md.
  -DF2DOCK_ENABLE_GPL_COMPONENTS=OFF
)

# The deps prefix carries libcvc (find_package(cvc CONFIG)), FFTW and the
# whole transitive closure libcvc's config file find_dependency()s.
if [[ -n "${CVC_DEPS_PREFIX:-}" ]]; then
  CMAKE_ARGS+=("-DCMAKE_PREFIX_PATH=$CVC_DEPS_PREFIX")
fi

cmake "${CMAKE_ARGS[@]}"
cmake --build "$CVC_BUILD_DIR" -j "$CVC_JOBS"
cmake --install "$CVC_BUILD_DIR"

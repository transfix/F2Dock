cmake_minimum_required(VERSION 3.16)
project(F2Dock VERSION 2.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Global Variables
if(MSVC)
  include_directories(${CMAKE_SOURCE_DIR}/inc)
  add_compile_options(/W3 /EHsc)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS -DNOMINMAX -D_USE_MATH_DEFINES)
  # POSIX strcasecmp is _stricmp on MSVC
  add_definitions(-Dstrcasecmp=_stricmp -Dstrncasecmp=_strnicmp)
else()
  set(CMAKE_CXX_FLAGS "$ENV{CPPFLAGS} -I${CMAKE_SOURCE_DIR}/inc -O2")
  set(CMAKE_LD_FLAGS "$ENV{LDFLAGS}")
endif()

find_package(Threads REQUIRED)
if(WIN32)
  find_package(PThreads4W REQUIRED)
endif()

# CMake Modules
SET(CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/CMake
  ${CMAKE_MODULE_PATH})

# include macro for setting up FFTW
include( SetupFFTW )

SetupFFTW()

# ── Optional GPL components ──────────────────────────────────────────
# F2Dock itself is LGPL-2.1-only (see LICENSE). Two vendored modules are
# not:
#
#   inc/fft-utils/sparsefft3*.h, src/fft-utils/sparsefft3*.cpp
#       GPL-2.0-or-later (MIT 2000, modified at CVC Lab) -- the sparse
#       3-D FFT used by the `useSparseFFT` parameter.
#   inc/libicp/*.h, src/f3dock-icp/{icp,icpPointToPoint,icpPointToPlane,
#   matrix}.cpp
#       GPL-2.0-or-later (Andreas Geiger, KIT) -- the ICP solver behind
#       the `icpRefine` parameter. src/f3dock-icp/kdtree.cpp is AFL-1.1
#       and references a LICENSE file that is not in the tree.
#
# Linking either into F2Dock makes the resulting binary a combined work
# covered by the GPL, which is not what this project distributes. They are
# therefore OFF by default. With the option OFF the sources are not
# compiled and not linked, `useSparseFFT` is forced off in favour of the
# dense FFTW path, and `icpRefine` is a no-op.
#
# Turning this ON produces a GPL-2.0-or-later binary. Do not ship the
# result under F2Dock's own licence. See THIRD-PARTY.md.
option(F2DOCK_ENABLE_GPL_COMPONENTS
  "Build the GPL-licensed sparse FFT and libicp modules. The resulting binary is covered by the GPL, not the LGPL."
  OFF)

if(F2DOCK_ENABLE_GPL_COMPONENTS)
  message(WARNING
    "F2DOCK_ENABLE_GPL_COMPONENTS=ON: sparsefft3 and libicp are GPL-2.0-or-later. "
    "The F2Dock binary produced by this build is a combined work covered by the "
    "GNU GPL v2 or later, not by F2Dock's own LGPL-2.1. See THIRD-PARTY.md.")
else()
  message(STATUS
    "F2Dock: GPL components disabled (default). Sparse FFT and ICP refinement "
    "are unavailable; the dense FFTW transform is used. "
    "Re-run with -DF2DOCK_ENABLE_GPL_COMPONENTS=ON to build them.")
endif()

option(F2DOCK_ENABLE_F3DOCK_EXPERIMENTAL "Build experimental features ported from F3Dock" ON)
option(F2DOCK_ENABLE_F3DOCK_SCIENCE_STACK "Discover optional scientific dependencies for F3Dock ports" ON)

if(F2DOCK_ENABLE_F3DOCK_EXPERIMENTAL)
  add_library(f3dock_science_stack INTERFACE)
  if(F2DOCK_ENABLE_F3DOCK_SCIENCE_STACK)
    find_package(Eigen3 QUIET)
    if(Eigen3_FOUND)
      target_link_libraries(f3dock_science_stack INTERFACE Eigen3::Eigen)
      target_compile_definitions(f3dock_science_stack INTERFACE F2DOCK_HAVE_EIGEN3=1)
    endif()

    find_package(LAPACK QUIET)
    if(LAPACK_FOUND)
      if(TARGET LAPACK::LAPACK)
        target_link_libraries(f3dock_science_stack INTERFACE LAPACK::LAPACK)
      elseif(LAPACK_LIBRARIES)
        target_link_libraries(f3dock_science_stack INTERFACE ${LAPACK_LIBRARIES})
      endif()
      target_compile_definitions(f3dock_science_stack INTERFACE F2DOCK_HAVE_LAPACK=1)
    endif()

    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
      pkg_check_modules(NFFT QUIET IMPORTED_TARGET nfft3)
      if(NFFT_FOUND)
        target_link_libraries(f3dock_science_stack INTERFACE PkgConfig::NFFT)
        target_compile_definitions(f3dock_science_stack INTERFACE F2DOCK_HAVE_NFFT=1)
        message(STATUS "F2Dock: NFFT3 ${NFFT_VERSION} detected via pkg-config; F2DOCK_HAVE_NFFT enabled")
      else()
        message(STATUS "F2Dock: NFFT3 not found via pkg-config; F2DOCK_HAVE_NFFT disabled")
      endif()
    endif()
  endif()
endif()

# ── libcvc dependency ────────────────────────────────────────────────
# libcvc comes from cvcpkg. Nothing else: no FetchContent source build, no
# release-tarball download, no system package manager.
#
# Provision the prefix from F2Dock's own recipe, which declares the closure:
#
#     cvcpkg install-deps cvcpkg/recipes/f2dock --prefix deps \
#         --config release --link shared
#     cmake -B build -DCMAKE_PREFIX_PATH="$PWD/deps"
#
# That installs libcvc plus everything cvcConfig.cmake find_dependency()s
# (Boost, assimp) and everything F2Dock links directly (FFTW3, pthreads4w on
# Windows) -- ~40 components, all hermetic cvcpkg bundles.
#
# Why the source build is gone rather than kept as a fallback: compiling
# libcvc here is what pulled CGAL and the host toolchain into F2Dock's build,
# and it broke on both non-Linux platforms at once. CGAL 6.2 (which Homebrew
# and vcpkg both ship now) collides with a `#define error(x)` in libcvc's
# reg3data.h, and windows-latest's VS18 v180 crashes CL.exe outright compiling
# cvc.vcxproj. A published bundle compiles none of that. A fallback that only
# runs when the fast path fails is a fallback that only runs on the broken
# platforms, so there is no version of it worth keeping.
set(LIBCVC_VERSION "3.3.0" CACHE STRING "Minimum libcvc version to require")

find_package(cvc ${LIBCVC_VERSION} CONFIG QUIET)

if(NOT cvc_FOUND)
  message(FATAL_ERROR
    "libcvc ${LIBCVC_VERSION} or newer was not found.\n"
    "F2Dock consumes libcvc from cvcpkg. Provision it with:\n"
    "    cvcpkg install-deps cvcpkg/recipes/f2dock --prefix deps "
    "--config release --link shared\n"
    "then configure with -DCMAKE_PREFIX_PATH=<abs path to deps>.")
endif()

# F2DockServer is built on libcvc's XML-RPC module. A libcvc bundle published
# without -DCVC_USING_XMLRPC=ON exports no cvc::xmlrpc and ships no
# lib/libxmlrpc archive, only include/xmlrpc/ headers -- fail loudly and name
# the minimum revision rather than letting the link fail later.
if(NOT TARGET cvc::xmlrpc)
  message(FATAL_ERROR
    "libcvc ${cvc_VERSION} was found but exports no cvc::xmlrpc target.\n"
    "F2DockServer requires it. Bundles before libcvc 3.3.0+cvc.4 were "
    "published with CVC_USING_XMLRPC off and are unusable here; install "
    "3.3.0+cvc.4 or newer.")
endif()

message(STATUS "Using libcvc ${cvc_VERSION} from ${cvc_DIR}")
#output 
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
mark_as_advanced(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

IF(APPLE)
ADD_DEFINITIONS(-DMACOS_X)
ENDIF(APPLE)

# ── Version header ──
configure_file(
  "${CMAKE_SOURCE_DIR}/inc/f2dock/Version.h.in"
  "${CMAKE_BINARY_DIR}/inc/f2dock/Version.h"
  @ONLY)
include_directories("${CMAKE_BINARY_DIR}/inc")

# Applied here rather than beside the option() above so the definitions reach
# F2Dock's own targets (src/ and tests/unit/, which compiles IcpRefine.cpp
# directly) without leaking into the libcvc subproject that FetchContent adds
# between the two points.
if(F2DOCK_ENABLE_GPL_COMPONENTS)
  add_compile_definitions(F2DOCK_WITH_SPARSEFFT=1 F2DOCK_WITH_LIBICP=1)
endif()

add_subdirectory (src)

# ── Testing ──
option(F2DOCK_BUILD_TESTS "Build unit tests" ON)
if(F2DOCK_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests/unit)
endif()

option(F2DOCK_BUILD_INTEGRATION_TESTS "Build long-running integration tests" OFF)
if(F2DOCK_BUILD_TESTS AND F2DOCK_BUILD_INTEGRATION_TESTS)
  add_subdirectory(tests/integration)
endif()

# ── Install targets ──
set(_F2DOCK_INSTALL_TARGETS F2Dock GB-Rerank FilterOutput)
if(TARGET F2DockServer)
  list(APPEND _F2DOCK_INSTALL_TARGETS F2DockServer)
endif()
install(TARGETS ${_F2DOCK_INSTALL_TARGETS}
        RUNTIME DESTINATION bin)

# Install data files needed at runtime
install(FILES
  deg04.matrix deg06.matrix deg08.matrix deg10.matrix
  deg12.matrix deg20.matrix
  DESTINATION share/F2Dock/data)

# Install documentation and examples
install(FILES README.md LICENSE DESTINATION .)
install(FILES
  tests/1A2K.inp tests/1ACB.inp tests/1AHW.inp tests/1JZD.inp
  DESTINATION share/F2Dock/examples)
install(DIRECTORY tests/
  DESTINATION share/F2Dock/examples
  FILES_MATCHING
    PATTERN "*.f2d"
    PATTERN "*.pqr"
    PATTERN "*.quad"
    PATTERN "*.txt"
    PATTERN "unit" EXCLUDE)

# ── CPack configuration ──
set(CPACK_PACKAGE_NAME "F2Dock")
set(CPACK_PACKAGE_VENDOR "CVC, The University of Texas at Austin")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Fast Fourier protein-protein docking")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")

# DEB-specific
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "CVC Lab <bajaj@cs.utexas.edu>")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libfftw3-3")

include(CPack)

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

include(FetchContent)
# ── libcvc dependency ────────────────────────────────────────────────
# Resolution order, falling through on failure:
#   1. find_package(cvc CONFIG) against a pre-installed libcvc on
#      CMAKE_PREFIX_PATH (fastest, no download or build).
#   2. Download the prebuilt release tarball matching this platform /
#      build type and find_package against the extracted tree.
#   3. FetchContent_Declare(GIT_REPOSITORY ...) and build from source.
#
# F2Dock links cvc::xmlrpc (XML-RPC support used by F2DockServer). When
# the prebuilt archive does not ship cvc::xmlrpc, or when no prebuilt
# archive has been published for the requested tag at all, strategies
# 1 and 2 are rejected and we fall through to the source build (which
# in turn pulls libcvc-deps for the heavy science stack, including
# NFFT). The source path is the supported way to consume a libcvc
# release that does not yet have platform binaries attached.
set(LIBCVC_TAG     "v3.2.1" CACHE STRING "libcvc release tag to consume")
# Package version that libcvc reports via cvcConfigVersion.cmake. The v3.2.1
# release artifacts still carry the 3.2.0 project version baked at build
# time (the bump shipped in #64 was tagged but not re-versioned in the
# CMake project() call). Keep the find_package() version match aligned
# with what's actually inside the prebuilt archive so the prebuilt path
# does not silently fall through to a from-source build.
set(LIBCVC_VERSION "3.2.0"  CACHE STRING "libcvc package version to look up")
set(LIBCVC_ARCHIVE_VERSION "${LIBCVC_VERSION}"
    CACHE STRING "Version string embedded in the prebuilt archive filename")
set(LIBCVC_GIT_SHA "08c15b5e244564469a34cddb6cdea00c87add906"
    CACHE STRING "libcvc commit SHA used when falling through to source build")
option(LIBCVC_FORCE_SOURCE "Skip prebuilt libcvc and always build from source" OFF)
set(LIBCVC_NEEDS_XMLRPC TRUE)

# Lean on libcvc for XML-RPC support used by F2DockServer.
set(CVC_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(CVC_ENABLE_CUDA OFF CACHE BOOL "" FORCE)
set(CVC_ENABLE_COVERAGE OFF CACHE BOOL "" FORCE)
set(CVC_BUILD_VOLROVER3 OFF CACHE BOOL "" FORCE)
set(CVC_USING_XMLRPC ON CACHE BOOL "" FORCE)
# F2Dock does not load MRC volumes. Disable libcvc's IMOD MRC path so the
# source build does not require libiimod (only shipped via libcvc-deps as
# of libcvc 3.2.1).
set(CVC_USING_IMOD_MRC OFF CACHE BOOL "" FORCE)

function(_libcvc_check_usable out_var)
  set(${out_var} FALSE PARENT_SCOPE)
  if(NOT cvc_FOUND)
    return()
  endif()
  if(LIBCVC_NEEDS_XMLRPC AND NOT TARGET cvc::xmlrpc)
    message(STATUS "libcvc package found (${cvc_VERSION}) but cvc::xmlrpc is missing")
    return()
  endif()
  set(${out_var} TRUE PARENT_SCOPE)
endfunction()

set(_libcvc_usable FALSE)

# 1. Pre-installed package
if(NOT LIBCVC_FORCE_SOURCE)
  find_package(cvc ${LIBCVC_VERSION} CONFIG QUIET)
  _libcvc_check_usable(_libcvc_usable)
  if(_libcvc_usable)
    message(STATUS "Using pre-installed libcvc ${cvc_VERSION}")
  endif()
endif()

# 2. Prebuilt release tarball
if(NOT _libcvc_usable AND NOT LIBCVC_FORCE_SOURCE)
  if(WIN32)
    set(_libcvc_plat "windows-x86_64")
    set(_libcvc_ext  "zip")
  elseif(APPLE)
    set(_libcvc_plat "macos-arm64")
    set(_libcvc_ext  "zip")
  elseif(UNIX)
    set(_libcvc_plat "linux-x86_64")
    set(_libcvc_ext  "tar.gz")
  endif()
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(_libcvc_cfg "debug")
  else()
    set(_libcvc_cfg "release")
  endif()
  if(DEFINED _libcvc_plat)
    set(_libcvc_archive "libcvc-${LIBCVC_ARCHIVE_VERSION}-${_libcvc_plat}-${_libcvc_cfg}.${_libcvc_ext}")
    set(_libcvc_url     "https://github.com/transfix/libcvc/releases/download/${LIBCVC_TAG}/${_libcvc_archive}")
    set(_libcvc_stage   "${CMAKE_BINARY_DIR}/libcvc-prebuilt")
    set(_libcvc_cfgfile "${_libcvc_stage}/lib/cmake/cvc/cvcConfig.cmake")
    if(NOT EXISTS "${_libcvc_cfgfile}")
      message(STATUS "Downloading prebuilt libcvc: ${_libcvc_url}")
      file(MAKE_DIRECTORY "${_libcvc_stage}")
      file(DOWNLOAD "${_libcvc_url}" "${_libcvc_stage}/${_libcvc_archive}"
           STATUS _libcvc_dl TLS_VERIFY ON)
      list(GET _libcvc_dl 0 _libcvc_dl_code)
      if(_libcvc_dl_code EQUAL 0)
        execute_process(
          COMMAND ${CMAKE_COMMAND} -E tar xzf "${_libcvc_archive}"
          WORKING_DIRECTORY "${_libcvc_stage}"
          RESULT_VARIABLE _libcvc_extract)
        if(NOT _libcvc_extract EQUAL 0)
          message(STATUS "Failed to extract ${_libcvc_archive} (exit ${_libcvc_extract})")
        endif()
      else()
        message(STATUS "Prebuilt libcvc download failed: ${_libcvc_dl}")
      endif()
    endif()
    if(EXISTS "${_libcvc_cfgfile}")
      list(PREPEND CMAKE_PREFIX_PATH "${_libcvc_stage}")
      find_package(cvc ${LIBCVC_VERSION} CONFIG QUIET)
      _libcvc_check_usable(_libcvc_usable)
      if(_libcvc_usable)
        message(STATUS "Using prebuilt libcvc ${cvc_VERSION} from ${_libcvc_url}")
      endif()
    endif()
  endif()
endif()

# 3. Source build (fallback)
if(NOT _libcvc_usable)
  message(STATUS "Building libcvc ${LIBCVC_TAG} from source")
  FetchContent_Declare(
    libcvc
    GIT_REPOSITORY https://github.com/transfix/libcvc.git
    GIT_TAG        ${LIBCVC_GIT_SHA}
  )
  FetchContent_MakeAvailable(libcvc)
endif()

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
install(FILES README.md DESTINATION .)
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
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/INSTALL")

# DEB-specific
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "CVC Lab <bajaj@cs.utexas.edu>")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libfftw3-3")

include(CPack)

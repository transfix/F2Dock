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
      target_link_libraries(f3dock_science_stack INTERFACE LAPACK::LAPACK)
      target_compile_definitions(f3dock_science_stack INTERFACE F2DOCK_HAVE_LAPACK=1)
    endif()

    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
      pkg_check_modules(NFFT QUIET IMPORTED_TARGET nfft3)
      if(NFFT_FOUND)
        target_link_libraries(f3dock_science_stack INTERFACE PkgConfig::NFFT)
        target_compile_definitions(f3dock_science_stack INTERFACE F2DOCK_HAVE_NFFT=1)
      endif()
    endif()
  endif()
endif()

if(NOT WIN32)
  include(FetchContent)
  # Lean on libcvc for XML-RPC support used by F2DockServer.
  set(CVC_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(CVC_ENABLE_CUDA OFF CACHE BOOL "" FORCE)
  set(CVC_ENABLE_COVERAGE OFF CACHE BOOL "" FORCE)
  set(CVC_BUILD_VOLROVER3 OFF CACHE BOOL "" FORCE)
  set(CVC_USING_XMLRPC ON CACHE BOOL "" FORCE)
  FetchContent_Declare(
    libcvc
    GIT_REPOSITORY https://github.com/transfix/libcvc.git
    GIT_TAG        972afd8552b4755356e572f9b4ffc33d4366cda9 # v3.1.0
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
if(NOT WIN32)
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

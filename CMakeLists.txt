cmake_minimum_required(VERSION 3.16)
project(F2Dock LANGUAGES C CXX)

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

# CMake Modules
SET(CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/CMake
  ${CMAKE_MODULE_PATH})

# include macro for setting up FFTW
include( SetupFFTW )

SetupFFTW()

#output 
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
mark_as_advanced(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

IF(APPLE)
ADD_DEFINITIONS(-DMACOS_X)
ENDIF(APPLE)

add_subdirectory (src)

# ── Testing ──
option(F2DOCK_BUILD_TESTS "Build unit tests" ON)
if(F2DOCK_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests/unit)
endif()

# ── Install targets ──
install(TARGETS F2Dock GB-Rerank FilterOutput F2DockServer
        RUNTIME DESTINATION bin)

# Install data files needed at runtime
install(FILES
  deg04.matrix deg06.matrix deg08.matrix deg10.matrix
  deg12.matrix deg20.matrix
  DESTINATION share/F2Dock/data)

# ── CPack configuration ──
set(CPACK_PACKAGE_NAME "F2Dock")
set(CPACK_PACKAGE_VENDOR "CVC, The University of Texas at Austin")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Fast Fourier protein-protein docking")
set(CPACK_PACKAGE_VERSION_MAJOR 2)
set(CPACK_PACKAGE_VERSION_MINOR 0)
set(CPACK_PACKAGE_VERSION_PATCH 0)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/INSTALL")

# DEB-specific
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "CVC Lab <bajaj@cs.utexas.edu>")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libfftw3-3")

include(CPack)

cmake_minimum_required(VERSION 3.16)


FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../inc/f2dock/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

INCLUDE_DIRECTORIES(
  ../../inc
  ../../inc/f2dock
)


LINK_DIRECTORIES(
  $ENV{TACC_FFTW3_LIB}
)


SET(F2DOCK_LINK_LIBS 
  fast-clash
  fast-GB 
  fast-hydro
  fast-LJ 
  fast-PQ
  fast-resCont
  PG-range
  math
  utils
  vol
  fft-utils
  misc-ident
  # No ${FFTW_LIB} here: the transform library is whatever
  # F2DOCK_FFT_PROVIDER selected, and fft-utils links it PUBLIC so it
  # propagates. Naming FFTW here would drag the GPL library into every
  # build regardless of the provider chosen.
  Threads::Threads
)
if(UNIX)
  list(APPEND F2DOCK_LINK_LIBS c m)
endif()
if(WIN32)
  list(APPEND F2DOCK_LINK_LIBS PThreads4W::PThreads4W)
endif()

if(HAVE_LIBMOL)
SET(F2DOCK_LINK_LIBS  ${F2DOCK_LINK_LIBS} 
  hbondFilter
  libmol
)
endif(HAVE_LIBMOL)

if(TARGET f3dock-flex)
  list(APPEND F2DOCK_LINK_LIBS f3dock-flex)
endif()

if(TARGET f3dock-domain)
  list(APPEND F2DOCK_LINK_LIBS f3dock-domain)
endif()

if(TARGET f3dock-icp)
  list(APPEND F2DOCK_LINK_LIBS f3dock-icp)
endif()

# Phase 5 task 1: pull in the optional scientific stack INTERFACE
# target so F2Dock's startup diagnostic can compile-time detect
# F2DOCK_HAVE_NFFT (and future F2DOCK_HAVE_* macros). The INTERFACE
# target adds nothing to the link line when the optional deps are
# not present.
if(TARGET f3dock_science_stack)
  list(APPEND F2DOCK_LINK_LIBS f3dock_science_stack)
endif()


ADD_EXECUTABLE (F2Dock ${SOURCE_FILES} ${INCLUDE_FILES})


TARGET_LINK_LIBRARIES(F2Dock ${F2DOCK_LINK_LIBS})

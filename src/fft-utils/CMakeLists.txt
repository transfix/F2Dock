cmake_minimum_required(VERSION 3.16)

# sparsefft3.cpp / sparsefft3-plan.cpp are GPL-2.0-or-later; the rest of
# this directory is LGPL-2.1. Keep them out of the glob so the default build
# links no GPL object, and swap in the LGPL stub that satisfies the three
# entry points instead. See inc/fft-utils/sparsefft3-api.h.
SET(SOURCE_FILES
  fastfft.cpp
  rank-fftw.cpp
)
if(F2DOCK_ENABLE_GPL_COMPONENTS)
  LIST(APPEND SOURCE_FILES sparsefft3.cpp sparsefft3-plan.cpp)
else()
  LIST(APPEND SOURCE_FILES sparsefft3-stub.cpp)
endif()

FILE(GLOB INCLUDE_FILES ../../inc/fft-utils/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

INCLUDE_DIRECTORIES(
  ../../inc/fft-utils
  ../../inc
)

ADD_LIBRARY (fft-utils STATIC ${SOURCE_FILES} ${INCLUDE_FILES})

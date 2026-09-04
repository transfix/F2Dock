cmake_minimum_required(VERSION 3.16)

# sparsefft3.cpp / sparsefft3-plan.cpp are GPL-2.0-or-later; the rest of
# this directory is LGPL-2.1. Keep them out of the glob so the default build
# links no GPL object, and swap in the LGPL stub that satisfies the three
# entry points instead. See inc/fft-utils/sparsefft3-api.h.
SET(SOURCE_FILES
  fastfft.cpp
  # rank-fftw times candidate grid sizes to build fftw-rank.txt, and
  # computeEffGrid() is called from the engine, so it builds for every
  # provider -- which is also more correct: the efficient sizes should
  # reflect the transform actually in use. Its only FFTW-specific calls are
  # the wisdom pair, already dead behind #ifdef SAVE_WISDOM_TO_FILE (never
  # defined anywhere in the tree).
  rank-fftw.cpp
)

# Exactly one provider backend is compiled.
if(_f2dock_fft_provider STREQUAL "fftw")
  LIST(APPEND SOURCE_FILES fft_fftw.cpp)
elseif(_f2dock_fft_provider STREQUAL "pocketfft")
  LIST(APPEND SOURCE_FILES fft_pocketfft.cpp)
elseif(_f2dock_fft_provider STREQUAL "mkl")
  LIST(APPEND SOURCE_FILES fft_mkl.cpp)
elseif(_f2dock_fft_provider STREQUAL "cufft")
  LIST(APPEND SOURCE_FILES fft_cufft.cpp)
endif()
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

if(_f2dock_fft_provider STREQUAL "fftw")
  target_compile_definitions(fft-utils PUBLIC F2DOCK_FFT_PROVIDER_FFTW=1)
  target_link_libraries(fft-utils PUBLIC ${FFTW_LIB})
elseif(_f2dock_fft_provider STREQUAL "pocketfft")
  # Header-only, from the cvcpkg pocketfft bundle. find_path rather than
  # find_package: upstream ships no CMake config, only the header.
  find_path(POCKETFFT_INCLUDE_DIR pocketfft_hdronly.h REQUIRED)
  target_include_directories(fft-utils PUBLIC ${POCKETFFT_INCLUDE_DIR})
  message(STATUS "F2Dock: FFT provider = PocketFFT (BSD-3) from ${POCKETFFT_INCLUDE_DIR}")
elseif(_f2dock_fft_provider STREQUAL "mkl")
  find_package(MKL CONFIG REQUIRED)
  target_link_libraries(fft-utils PUBLIC MKL::MKL)
  message(STATUS "F2Dock: FFT provider = Intel oneMKL")
elseif(_f2dock_fft_provider STREQUAL "cufft")
  find_package(CUDAToolkit REQUIRED)
  target_link_libraries(fft-utils PUBLIC CUDA::cufft)
  message(STATUS "F2Dock: FFT provider = NVIDIA cuFFT")
endif()

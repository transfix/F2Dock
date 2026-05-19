cmake_minimum_required(VERSION 3.16)

# Phase 5 task 2: f3dock-nfft static library. Provides the FastGaussianSum
# wrapper that subsequent PRs will route the electrostatics gridding
# through under the `nfftAccelerate` opt-in. The library always builds (no
# external dependency on NFFT3 in the C++ surface); the NFFT fast path is
# selected at runtime via `FastGaussianSum::is_accelerated()` and is
# currently a no-op fallback to the direct O(N*M) evaluator.

add_library(f3dock-nfft STATIC
  FastGaussianSum.cpp
)

target_include_directories(f3dock-nfft
  PUBLIC
    ${CMAKE_SOURCE_DIR}/inc
)

# Link the optional INTERFACE target so F2DOCK_HAVE_NFFT (and any related
# include paths / link libs) propagates when libnfft3 was discovered at
# configure time. INTERFACE targets contribute nothing when the optional
# dependency is absent.
if(TARGET f3dock_science_stack)
  target_link_libraries(f3dock-nfft PUBLIC f3dock_science_stack)
endif()

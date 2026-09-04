cmake_minimum_required(VERSION 3.16)

# GoogleTest comes from cvcpkg like every other dependency -- it is declared
# in cvcpkg/recipes/f2dock/recipe.yaml, so `cvcpkg install-deps` puts it in
# the same prefix as libcvc and FFTW. This used to FetchContent v1.14.0 from
# github at configure time, which was the last network fetch in an otherwise
# hermetic build.
find_package(GTest CONFIG REQUIRED)

include(GoogleTest)

# Common include path for all test targets
set(F2DOCK_INC "${CMAKE_SOURCE_DIR}/inc")
set(F2DOCK_TESTS_DATA "${CMAKE_SOURCE_DIR}/tests")

# ── Helper: add a test executable that links against specific F2Dock libs ──
function(f2dock_add_test TEST_NAME)
  cmake_parse_arguments(ARG "" "" "SOURCES;LIBS" ${ARGN})
  add_executable(${TEST_NAME} ${ARG_SOURCES})
  target_include_directories(${TEST_NAME} PRIVATE ${F2DOCK_INC})
  target_link_libraries(${TEST_NAME} PRIVATE ${ARG_LIBS} GTest::gtest_main)
  target_compile_definitions(${TEST_NAME} PRIVATE
    F2DOCK_TEST_DATA_DIR="${F2DOCK_TESTS_DATA}"
  )
  gtest_discover_tests(${TEST_NAME})
endfunction()

# ── Math library tests (Matrix, Vector, Tuple) ──
f2dock_add_test(test_math
  SOURCES test_math.cpp
  LIBS math
)

# ── PG-range data structure tests ──
f2dock_add_test(test_pg_range
  SOURCES test_pg_range.cpp
  LIBS PG-range
)

# ── PairingHeap tests ──
f2dock_add_test(test_pairing_heap
  SOURCES test_pairing_heap.cpp
  LIBS fast-PQ
)

# ── clashFilter tests ──
set(CLASH_FILTER_LIBS fast-clash math Threads::Threads)
if(WIN32)
  list(APPEND CLASH_FILTER_LIBS PThreads4W::PThreads4W)
endif()
f2dock_add_test(test_clash_filter
  SOURCES test_clash_filter.cpp
  LIBS ${CLASH_FILTER_LIBS}
)

# ── File format tests (F2D, PQR, QUAD, parameter parsing) ──
f2dock_add_test(test_file_formats
  SOURCES test_file_formats.cpp
  LIBS utils
)

# ── ValuePosition3D tests ──
# ValuePosition3D.cpp lives in src/f2dock/ (it's part of the F2Dock executable,
# not a separate library), so we compile it directly into the test.
add_executable(test_value_position
  test_value_position.cpp
  ${CMAKE_SOURCE_DIR}/src/f2dock/ValuePosition3D.cpp
)
target_include_directories(test_value_position PRIVATE
  ${F2DOCK_INC}
  ${F2DOCK_INC}/f2dock   # ValuePosition3D.cpp includes "ValuePosition3D.h"
)
target_link_libraries(test_value_position PRIVATE math GTest::gtest_main)
target_compile_definitions(test_value_position PRIVATE
  F2DOCK_TEST_DATA_DIR="${F2DOCK_TESTS_DATA}"
)
gtest_discover_tests(test_value_position)

if(TARGET f3dock-loop-closure)
  f2dock_add_test(test_f3dock_loop_closure
    SOURCES test_f3dock_loop_closure.cpp
    LIBS f3dock-loop-closure
  )
endif()

if(TARGET f3dock-icp)
  f2dock_add_test(test_f3dock_icp
    SOURCES test_f3dock_icp.cpp
    LIBS f3dock-icp
  )
  f2dock_add_test(test_f2dock_icp_refine
    SOURCES
      test_f2dock_icp_refine.cpp
      ${CMAKE_SOURCE_DIR}/src/f2dock/IcpRefine.cpp
    LIBS f3dock-icp
  )
endif()

if(TARGET f3dock-flex)
  f2dock_add_test(test_f3dock_flex
    SOURCES test_f3dock_flex.cpp
    LIBS f3dock-flex
  )
  f2dock_add_test(test_f3dock_flex_sampler
    SOURCES test_f3dock_flex_sampler.cpp
    LIBS f3dock-flex
  )
  f2dock_add_test(test_f3dock_flex_apply
    SOURCES test_f3dock_flex_apply.cpp
    LIBS f3dock-flex
  )
endif()

if(TARGET f3dock-domain)
  f2dock_add_test(test_f3dock_domain
    SOURCES test_f3dock_domain.cpp
    LIBS f3dock-domain
  )
  f2dock_add_test(test_f3dock_domain_spec
    SOURCES test_f3dock_domain_spec.cpp
    LIBS f3dock-domain
  )
  f2dock_add_test(test_f3dock_domain_partition
    SOURCES test_f3dock_domain_partition.cpp
    LIBS f3dock-domain
  )
  f2dock_add_test(test_f3dock_domain_sampler
    SOURCES test_f3dock_domain_sampler.cpp
    LIBS f3dock-domain
  )
endif()

if(TARGET f3dock-nfft)
  f2dock_add_test(test_f3dock_nfft_fast_gaussian_sum
    SOURCES test_f3dock_nfft_fast_gaussian_sum.cpp
    LIBS f3dock-nfft
  )
endif()

# ── Header-only DockMode parser tests (no library dep) ──
add_executable(test_dock_mode test_dock_mode.cpp)
target_include_directories(test_dock_mode PRIVATE ${F2DOCK_INC})
target_link_libraries(test_dock_mode PRIVATE GTest::gtest_main)
gtest_discover_tests(test_dock_mode)

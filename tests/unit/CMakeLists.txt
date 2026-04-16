cmake_minimum_required(VERSION 3.16)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        f8d7d77c06936315286eb55f8de22cd23c188571  # v1.14.0
)
# Prevent GTest from overriding our compiler/linker options on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
# Don't install GTest headers/libs with our project
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

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

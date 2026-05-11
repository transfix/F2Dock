cmake_minimum_required(VERSION 3.16)

find_package(Python3 COMPONENTS Interpreter REQUIRED)

# This test is intentionally optional because it runs a full docking workflow
# and can take significantly longer than unit tests.
function(add_f2dock_integration_case case_name)
  cmake_parse_arguments(ARG "ALLOW_GENERATE_BASELINE" "INPUT;EXPECTED;ROWS;ABS_TOL;REL_TOL" "" ${ARGN})

  if(NOT EXISTS "${ARG_INPUT}")
    message(FATAL_ERROR "Missing integration input file: ${ARG_INPUT}")
  endif()
  if(ARG_ALLOW_GENERATE_BASELINE)
    set(_allow_generate ON)
  else()
    if(NOT EXISTS "${ARG_EXPECTED}")
      message(FATAL_ERROR "Missing integration expected output: ${ARG_EXPECTED}")
    endif()
    set(_allow_generate OFF)
  endif()

  add_test(
    NAME integration_${case_name}
    COMMAND ${CMAKE_COMMAND}
      -DF2DOCK_EXE=$<TARGET_FILE:F2Dock>
      -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
      -DBINARY_DIR=${CMAKE_BINARY_DIR}
      -DCASE_NAME=${case_name}
      -DINPUT_TEMPLATE=${ARG_INPUT}
      -DEXPECTED_OUTPUT=${ARG_EXPECTED}
      -DCOMPARE_ROWS=${ARG_ROWS}
      -DABS_TOL=${ARG_ABS_TOL}
      -DREL_TOL=${ARG_REL_TOL}
      -DALLOW_GENERATE_BASELINE=${_allow_generate}
      -DPYTHON_EXECUTABLE=${Python3_EXECUTABLE}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/RunAndCompareF2Dock.cmake
  )

  set_tests_properties(integration_${case_name} PROPERTIES
    LABELS "integration"
    TIMEOUT 1800
  )
endfunction()

add_f2dock_integration_case(
  1ACB
  INPUT    "${CMAKE_CURRENT_SOURCE_DIR}/inputs/1ACB_smoke.inp"
  EXPECTED "${CMAKE_CURRENT_SOURCE_DIR}/expected/1ACB_smoke_out.txt"
  ROWS     50
  ABS_TOL  1e-3
  REL_TOL  1e-5
)

add_f2dock_integration_case(
  1A2K
  INPUT    "${CMAKE_CURRENT_SOURCE_DIR}/inputs/1A2K_smoke.inp"
  EXPECTED "${CMAKE_CURRENT_SOURCE_DIR}/expected/1A2K_smoke_out.txt"
  ROWS     50
  ABS_TOL  1e-3
  REL_TOL  1e-5
)

if(TARGET F2DockServer)
  add_test(
    NAME integration_F2DockServer
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_f2dock_server.py
      --server-exe $<TARGET_FILE:F2DockServer>
      --work-dir ${CMAKE_BINARY_DIR}/tests/integration/f2dockserver_smoke
      --platform 0
  )

  set_tests_properties(integration_F2DockServer PROPERTIES
    LABELS "integration"
    TIMEOUT 120
  )
endif()

cmake_minimum_required(VERSION 3.16)

foreach(var F2DOCK_EXE SOURCE_DIR BINARY_DIR CASE_NAME INPUT_TEMPLATE EXPECTED_OUTPUT COMPARE_ROWS ABS_TOL REL_TOL ALLOW_GENERATE_BASELINE PYTHON_EXECUTABLE)
  if(NOT DEFINED ${var})
    message(FATAL_ERROR "Missing required variable: ${var}")
  endif()
endforeach()

if(NOT EXISTS "${F2DOCK_EXE}")
  message(FATAL_ERROR "F2Dock executable not found: ${F2DOCK_EXE}")
endif()
if(NOT EXISTS "${INPUT_TEMPLATE}")
  message(FATAL_ERROR "Input template not found: ${INPUT_TEMPLATE}")
endif()
if(NOT EXISTS "${PYTHON_EXECUTABLE}")
  message(FATAL_ERROR "Python executable not found: ${PYTHON_EXECUTABLE}")
endif()

set(work_dir "${BINARY_DIR}/tests/integration")
file(MAKE_DIRECTORY "${work_dir}")

set(actual_output "${work_dir}/${CASE_NAME}_out.actual.txt")
set(temp_input "${work_dir}/${CASE_NAME}.integration.inp")
set(generated_baseline "${work_dir}/${CASE_NAME}_out.generated_baseline.txt")
set(compare_script "${SOURCE_DIR}/tests/integration/compare_f2dock_outputs.py")

if(NOT EXISTS "${compare_script}")
  message(FATAL_ERROR "Comparison script not found: ${compare_script}")
endif()

file(READ "${INPUT_TEMPLATE}" input_text)
string(REGEX REPLACE "outFile[ \t]+[^\r\n]+" "outFile ${actual_output}" input_text "${input_text}")
file(WRITE "${temp_input}" "${input_text}")

execute_process(
  COMMAND "${F2DOCK_EXE}" "${temp_input}"
  WORKING_DIRECTORY "${SOURCE_DIR}"
  RESULT_VARIABLE run_rc
)

if(NOT run_rc EQUAL 0)
  message(FATAL_ERROR "F2Dock execution failed with code ${run_rc}")
endif()

if(NOT EXISTS "${actual_output}")
  message(FATAL_ERROR "F2Dock did not produce output file: ${actual_output}")
endif()

set(reference_output "${EXPECTED_OUTPUT}")
set(generate_requested FALSE)
if(ALLOW_GENERATE_BASELINE)
  if(ALLOW_GENERATE_BASELINE STREQUAL "ON" OR ALLOW_GENERATE_BASELINE STREQUAL "TRUE" OR ALLOW_GENERATE_BASELINE STREQUAL "1")
    set(generate_requested TRUE)
  endif()
endif()

set(use_generated_baseline FALSE)
if(generate_requested)
  set(source_expected_size 0)
  if(EXISTS "${EXPECTED_OUTPUT}")
    file(SIZE "${EXPECTED_OUTPUT}" source_expected_size)
  endif()
  if(source_expected_size GREATER 0)
    set(reference_output "${EXPECTED_OUTPUT}")
  else()
    file(COPY_FILE "${actual_output}" "${generated_baseline}" ONLY_IF_DIFFERENT)
    set(reference_output "${generated_baseline}")
    set(use_generated_baseline TRUE)
  endif()
else()
  if(NOT EXISTS "${EXPECTED_OUTPUT}")
    message(FATAL_ERROR "Expected output not found: ${EXPECTED_OUTPUT}")
  endif()
endif()

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" "${compare_script}"
    --expected "${reference_output}"
    --actual "${actual_output}"
    --rows "${COMPARE_ROWS}"
    --abs-tol "${ABS_TOL}"
    --rel-tol "${REL_TOL}"
  RESULT_VARIABLE compare_rc
  OUTPUT_VARIABLE compare_out
  ERROR_VARIABLE compare_err
)

if(NOT compare_rc EQUAL 0)
  message(FATAL_ERROR
    "Integration regression comparison failed for ${CASE_NAME}.\n"
    "STDOUT:\n${compare_out}\n"
    "STDERR:\n${compare_err}")
endif()

if(use_generated_baseline)
  message(WARNING
    "Expected baseline for ${CASE_NAME} is missing/empty; generated temporary baseline at ${generated_baseline}")
endif()

message(STATUS
  "Integration regression check passed for ${CASE_NAME}: ${compare_out}")

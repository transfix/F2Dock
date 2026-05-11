cmake_minimum_required(VERSION 3.16)

add_library(f3dock-loop-closure STATIC
  tripep_closure.cpp
  LoopClosureSolver.cpp
)

target_include_directories(f3dock-loop-closure
  PUBLIC
    ${CMAKE_SOURCE_DIR}/inc
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

if(TARGET f3dock_science_stack)
  target_link_libraries(f3dock-loop-closure PUBLIC f3dock_science_stack)
endif()

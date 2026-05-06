cmake_minimum_required(VERSION 3.16)

add_library(f3dock-icp STATIC
  IcpAligner.cpp
  matrix.cpp
  kdtree.cpp
  icp.cpp
  icpPointToPoint.cpp
)

target_include_directories(f3dock-icp
  PUBLIC
    ${CMAKE_SOURCE_DIR}/inc
)

if(TARGET f3dock_science_stack)
  target_link_libraries(f3dock-icp PUBLIC f3dock_science_stack)
endif()

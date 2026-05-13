cmake_minimum_required(VERSION 3.16)

add_library(f3dock-flex STATIC
  FlexKinematics.cpp
  FlexSampler.cpp
)

target_include_directories(f3dock-flex
  PUBLIC
    ${CMAKE_SOURCE_DIR}/inc
)

if(TARGET f3dock_science_stack)
  target_link_libraries(f3dock-flex PUBLIC f3dock_science_stack)
endif()

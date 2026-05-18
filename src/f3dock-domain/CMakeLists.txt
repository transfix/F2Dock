cmake_minimum_required(VERSION 3.16)

add_library(f3dock-domain STATIC
  DomainGraph.cpp
  DomainPartition.cpp
  DomainSpec.cpp
)

target_include_directories(f3dock-domain
  PUBLIC
    ${CMAKE_SOURCE_DIR}/inc
)

if(TARGET f3dock-flex)
  target_link_libraries(f3dock-domain PUBLIC f3dock-flex)
endif()

if(TARGET f3dock_science_stack)
  target_link_libraries(f3dock-domain PUBLIC f3dock_science_stack)
endif()

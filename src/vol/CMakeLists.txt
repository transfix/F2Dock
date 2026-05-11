cmake_minimum_required(VERSION 3.16)

FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../inc/vol/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

INCLUDE_DIRECTORIES(
  ../../inc/vol
)

ADD_LIBRARY (vol STATIC ${SOURCE_FILES} ${INCLUDE_FILES})

# Volume serialization (RAWIV.cpp) now delegates to libcvc's volume I/O.
if(NOT TARGET cvc::cvc)
  message(FATAL_ERROR "cvc::cvc target not found; libcvc dependency was not configured")
endif()
target_link_libraries(vol PUBLIC cvc::cvc)

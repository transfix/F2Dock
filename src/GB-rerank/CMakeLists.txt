cmake_minimum_required(VERSION 3.16)


FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../inc/fGB-rerank/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

INCLUDE_DIRECTORIES(
  ../../inc
)

SET(GBRERANK_LINK_LIBS 
  fast-GB 
  PG-range
  utils
  Threads::Threads
)
if(UNIX)
  list(APPEND GBRERANK_LINK_LIBS c m)
endif()
if(WIN32)
  list(APPEND GBRERANK_LINK_LIBS PThreads4W::PThreads4W)
endif()


ADD_EXECUTABLE (GB-Rerank ${SOURCE_FILES} ${INCLUDE_FILES})


TARGET_LINK_LIBRARIES(GB-Rerank ${GBRERANK_LINK_LIBS})

cmake_minimum_required (VERSION 2.8.10 FATAL_ERROR)

message("CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message("INCLUDE_DIR: ${CMAKE_SOURCE_DIR}/../libs/libpll-2/src")

include_directories("${CMAKE_SOURCE_DIR}/../libs/libpll-2/src")

link_directories(${CMAKE_SOURCE_DIR}/../build/libs/libpll-2/src)
message("LINK_DIR: ${CMAKE_SOURCE_DIR}/../build/libs/libpll-2/src")
find_library(PLL_LIBRARY pll PATHS ${CMAKE_SOURCE_DIR}/../build/libs/libpll-2/src)
if(NOT PLL_LIBRARY)
    message(FATAL_ERROR "pll library not found")
endif()

add_executable(rtree_parser rtree_parser.c)
target_link_libraries(rtree_parser ${PLL_LIBRARY})

add_executable(init init.c init.h map.c)
target_link_libraries(init ${PLL_LIBRARY})

link_directories(${CMAKE_SOURCE_DIR}/../build/libs/libpll-2/src)
message("LINK_DIR: ${CMAKE_SOURCE_DIR}/../build/libs/libpll-2/src")
find_library(PLL_LIBRARY pll PATHS ${CMAKE_SOURCE_DIR}/../build/libs/libpll-2/src)
if(NOT PLL_LIBRARY)
    message(FATAL_ERROR "pll library not found")
endif()

target_link_libraries(rtree_parser ${PLL_LIBRARY})
target_link_libraries(init ${PLL_LIBRARY})

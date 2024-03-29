cmake_minimum_required (VERSION 2.8.10 FATAL_ERROR)

project (libqfml-debug)

if (NOT DEFINED PLLMODULES_LIBPLL_PATH)
    set(PLLMODULES_LIBPLL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/libs/libpll-2")
endif()

if(NOT DEFINED BUILD_PLLMODULES_SHARED AND NOT DEFINED BUILD_PLLMODULES_STATIC)
    set(BUILD_PLLMODULES_SHARED ${BUILD_SHARED_LIBS})
    if (BUILD_SHARED_LIBS)
        set(BUILD_PLLMODULES_STATIC FALSE)
    else()
        set(BUILD_PLLMODULES_STATIC TRUE)
    endif()
endif()


if(BUILD_PLLMODULES_SHARED)
    message(STATUS "pll-modules shared build enabled")
endif()
if(BUILD_PLLMODULES_STATIC)
    message(STATUS "pll-modules static build enabled")
endif()

if (NOT PLLMODULES_SKIP_LIBPLL_BUILD)
    set(BUILD_LIBPLL_SHARED ${BUILD_PLLMODULES_SHARED})
    set(BUILD_LIBPLL_STATIC ${BUILD_PLLMODULES_STATIC})
    add_subdirectory(${PLLMODULES_LIBPLL_PATH})
endif()

# understand what's going on here
macro(add_pllmodules_lib target sources)
    add_library(${target}_obj OBJECT ${sources})
    set_property(TARGET ${target}_obj PROPERTY POSITION_INDEPENDENT_CODE 1)
    target_include_directories(${target}_obj PRIVATE ${PLLMOD_INCLUDE_PATH})
    if(BUILD_PLLMODULES_STATIC)
        add_library(${target}_static STATIC $<TARGET_OBJECTS:${target}_obj>)
        set_target_properties(${target}_static PROPERTIES OUTPUT_NAME "${target}")
        target_link_libraries(${target}_static ${PLL_LIBRARIES})
        set(PLLMODULES_LIBRARIES
                ${target}_static ${PLLMODULES_LIBRARIES}
                CACHE INTERNAL "${PROJECT_NAME}: Libraries to link against")
    endif()
endmacro()

set(PLLMOD_EXPORTED_INCLUDE_PATH ${PLLMOD_EXPORTED_INCLUDE_PATH} ${PLL_EXPORTED_INCLUDE_PATH} PARENT_SCOPE)

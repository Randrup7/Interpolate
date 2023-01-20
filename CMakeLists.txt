cmake_minimum_required(VERSION 3.24.2)

project(Main 
        VERSION 1.0.0
        DESCRIPTION "optional description of project"
        LANGUAGES CXX )

set(CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

add_compile_options(
       -Wall -Werror 
       $<$<CONFIG:RELEASE>:-O2>
       $<$<CONFIG:DEBUG>:-O0>
       $<$<CONFIG:DEBUG>:-ggdb3>
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/exe")

add_subdirectory(src)
add_subdirectory(exe)
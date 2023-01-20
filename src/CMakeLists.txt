add_library(
    interpolateLib STATIC
    interpolate.cpp
)

# Specify library's header file location
target_include_directories(
    interpolateLib PUBLIC
    "${CMAKE_SOURCE_DIR}/include"
)

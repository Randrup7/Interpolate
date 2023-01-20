add_executable(${PROJECT_NAME} main.cpp)

# Linking library target (helloLib) to executable target ${PROJECT_NAME}
target_link_libraries(${PROJECT_NAME} PRIVATE interpolateLib)
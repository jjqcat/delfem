add_executable(scalar3d main.cpp)
link_directories("${PROJECT_SOURCE_DIR}/lib")
target_link_libraries(scalar3d delfemlib)
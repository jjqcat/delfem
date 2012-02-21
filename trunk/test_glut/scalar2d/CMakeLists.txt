add_executable(scalar2d main.cpp)
link_directories("${PROJECT_SOURCE_DIR}/lib")
target_link_libraries(scalar2d delfemlib)
add_executable(solid2d main.cpp)
link_directories("${PROJECT_SOURCE_DIR}/lib")
target_link_libraries(solid2d delfemlib)
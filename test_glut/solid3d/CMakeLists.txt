add_executable(solid3d main.cpp)
link_directories("${PROJECT_SOURCE_DIR}/lib")
target_link_libraries(solid3d delfemlib)
add_executable(fluid2d main.cpp)
link_directories("${PROJECT_SOURCE_DIR}/lib")
target_link_libraries(fluid2d delfemlib)
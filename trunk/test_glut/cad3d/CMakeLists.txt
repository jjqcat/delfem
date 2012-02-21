add_executable(cad3d main.cpp)
link_directories("${PROJECT_SOURCE_DIR}/lib")
target_link_libraries(cad3d delfemlib)
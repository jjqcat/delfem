add_executable(msh_view main.cpp)
link_directories("${PROJECT_SOURCE_DIR}/lib")
target_link_libraries(msh_view delfemlib)
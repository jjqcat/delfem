add_executable(cad2d_edit main.cpp)
link_directories("${PROJECT_SOURCE_DIR}/lib")
target_link_libraries(cad2d_edit delfemlib)
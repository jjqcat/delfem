add_executable(field_evaluate main.cpp)
link_directories("${PROJECT_SOURCE_DIR}/lib")
target_link_libraries(field_evaluate delfemlib)
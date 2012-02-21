cmake_minimum_required(VERSION 2.8)
project(DelFEM)
include_directories(include)
#add_subdirectory(lib  test_glut/cad2d)
subdirs(lib 
test_glut/cad2d 
test_glut/cad2d_edit
test_glut/cad3d
test_glut/msh_view
test_glut/field_evaluate
test_glut/scalar2d 
test_glut/scalar3d
test_glut/solid2d
test_glut/solid3d )
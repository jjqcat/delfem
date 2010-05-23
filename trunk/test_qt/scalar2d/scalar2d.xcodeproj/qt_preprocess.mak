#############################################################################
# Makefile for building: scalar2d.app/Contents/MacOS/scalar2d
# Generated by qmake (2.01a) (Qt 4.6.2) on: Sun May 23 10:25:56 2010
# Project:  scalar2d.pro
# Template: app
# Command: /usr/bin/qmake -spec /usr/local/Qt4.6/mkspecs/macx-xcode -macx -o scalar2d.xcodeproj/project.pbxproj scalar2d.pro
#############################################################################

MOC       = /Developer/Tools/Qt/moc
UIC       = /Developer/Tools/Qt/uic
LEX       = flex
LEXFLAGS  = 
YACC      = yacc
YACCFLAGS = -d
DEFINES       = -DQT_OPENGL_LIB -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED
INCPATH       = -I/usr/local/Qt4.6/mkspecs/macx-xcode -I. -I/Library/Frameworks/QtCore.framework/Versions/4/Headers -I/usr/include/QtCore -I/Library/Frameworks/QtGui.framework/Versions/4/Headers -I/usr/include/QtGui -I/Library/Frameworks/QtOpenGL.framework/Versions/4/Headers -I/usr/include/QtOpenGL -I/usr/include -I/System/Library/Frameworks/OpenGL.framework/Versions/A/Headers -I/System/Library/Frameworks/AGL.framework/Headers -I. -I. -I/usr/local/include -I/System/Library/Frameworks/CarbonCore.framework/Headers -I../../include -F/Library/Frameworks
DEL_FILE  = rm -f
MOVE      = mv -f

IMAGES = 
PARSERS =
preprocess: $(PARSERS) compilers
clean preprocess_clean: parser_clean compiler_clean

parser_clean:
mocclean: compiler_moc_header_clean compiler_moc_source_clean

mocables: compiler_moc_header_make_all compiler_moc_source_make_all

compilers: ./moc_mainwindow.cpp ./moc_glwidget.cpp ./ui_mainwindow.h
compiler_objective_c_make_all:
compiler_objective_c_clean:
compiler_moc_header_make_all: moc_mainwindow.cpp moc_glwidget.cpp
compiler_moc_header_clean:
	-$(DEL_FILE) moc_mainwindow.cpp moc_glwidget.cpp
moc_mainwindow.cpp: mainwindow.h
	/Developer/Tools/Qt/moc $(DEFINES) $(INCPATH) -D__APPLE__ -D__GNUC__ mainwindow.h -o moc_mainwindow.cpp

moc_glwidget.cpp: glwidget.h
	/Developer/Tools/Qt/moc $(DEFINES) $(INCPATH) -D__APPLE__ -D__GNUC__ glwidget.h -o moc_glwidget.cpp

compiler_rcc_make_all:
compiler_rcc_clean:
compiler_image_collection_make_all: qmake_image_collection.cpp
compiler_image_collection_clean:
	-$(DEL_FILE) qmake_image_collection.cpp
compiler_moc_source_make_all:
compiler_moc_source_clean:
compiler_rez_source_make_all:
compiler_rez_source_clean:
compiler_uic_make_all: ui_mainwindow.h
compiler_uic_clean:
	-$(DEL_FILE) ui_mainwindow.h
ui_mainwindow.h: mainwindow.ui
	/Developer/Tools/Qt/uic mainwindow.ui -o ui_mainwindow.h

compiler_yacc_decl_make_all:
compiler_yacc_decl_clean:
compiler_yacc_impl_make_all:
compiler_yacc_impl_clean:
compiler_lex_make_all:
compiler_lex_clean:
compiler_clean: compiler_moc_header_clean compiler_uic_clean 


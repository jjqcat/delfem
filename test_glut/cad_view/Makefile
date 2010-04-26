CXX    = g++
CFLAGS = -Wall -O2
LDFLAGS =
INCLUDES = -I../../include
LIBS = -L../../lib -ldfm
#Linux(defalut)
ifeq ($(OS),linux)
	LIBS_GL = -lglut -lGL -lGLU
else ifeq ($(OS),Windows_NT)
	LIBS_GL = -lfreeglut -lglu32 -lopengl32	# on Windows
else ifeq ($(shell uname -s),Darwin)
	LIBS_GL = -framework OpenGL -framework GLUT
endif
OBJS = main.o
TARGET = main.exe

all: $(TARGET)
					
$(TARGET): $(OBJS)
	$(OS)
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LIBS_GL) $(LIBS)

clean:
	-rm -f $(OBJS)
.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<
CXX    = g++
CFLAGS = -Wall -O2
LDFLAGS =
INCLUDES =
LIBS = -ldfm
ifeq ($(OS),Windows_NT)
	LIBS_GL = -lglut32 -lglu32 -lopengl32	# on Windows
else
	LIBS_GL = -lglut -lGL -lGLU
endif
OBJS = main.o
TARGET = main.exe

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LIBS) $(LIBS_GL)

clean:
	-rm -f $(OBJS)
.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<
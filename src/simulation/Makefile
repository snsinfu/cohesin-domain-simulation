SOURCES = $(wildcard *.cpp misc/*.cpp reactions/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

CXXFLAGS = \
  $(DBGFLAGS) \
  $(WARNFLAGS) \
  $(OPTFLAGS) \
  $(INCLUDES) \
  $(EXTRA_CXXFLAGS)

DBGFLAGS = \
  -g

WARNFLAGS = \
  -std=c++23 \
  -Wpedantic \
  -Wall \
  -Wextra \
  -Wconversion \
  -Wsign-conversion \
  -Wshadow

OPTFLAGS = \
  -O2 \
  -funsafe-math-optimizations

INCLUDES = \
  -isystem ../../submodules/github.com/snsinfu/micromd/include \
  -isystem ../../submodules/github.com/snsinfu/cxx-spline/include \
  -isystem ../../submodules/github.com/snsinfu/cxx-getopt \
  -isystem ../../submodules/github.com/danielaparker/jsoncons/include

TARGET = main

ARTIFACTS = $(TARGET) $(OBJECTS) _depends.mk


.PHONY: all clean depends

all: $(TARGET)
	@:

clean:
	@rm -f $(ARTIFACTS)

depends:
	for src in $(SOURCES); do \
	  $(CXX) $(CXXFLAGS) -c -MM -MF- -MT $${src%.*}.o $${src}; \
	done > _depends.mk

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

-include _depends.mk
-include _overrides.mk

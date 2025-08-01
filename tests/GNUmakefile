HEADERS := $(wildcard ../*.H)
SOURCES := $(wildcard *.cpp)
EXECUTABLES := $(SOURCES:.cpp=)

#DEFINES := -DUSE_LONG_DOUBLE
DEFINES := -DUSE_FLOAT128

INCLUDES := ..

INCLUDE_FLAGS := $(patsubst %,-I%,$(INCLUDES))
CXXFLAGS := -fext-numeric-literals -O3 -std=c++23 -Wall -Wextra -Wshadow

ALL: $(EXECUTABLES)

% : %.cpp $(HEADERS)
	g++  $(CXXFLAGS) $(DEFINES) $(INCLUDE_FLAGS) -o $@ $<

clean:
	rm -f $(EXECUTABLES)

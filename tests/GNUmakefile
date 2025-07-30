HEADERS := $(wildcard ../*.H)
SOURCES := $(wildcard *.cpp)
EXECUTABLES := $(SOURCES:.cpp=)

#DEFINES := -DUSE_LONG_DOUBLE
DEFINES := -DUSE_FLOAT128

ALL: $(EXECUTABLES)

% : %.cpp $(HEADERS)
	g++ -O3 -std=c++23 -Wall -Wextra -Wshadow -Wpedantic $(DEFINES) -I.. -o $@ $<

clean:
	rm -f $(EXECUTABLES)

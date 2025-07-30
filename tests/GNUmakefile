HEADERS := $(wildcard ../*.H)
SOURCES := $(wildcard *.cpp)
EXECUTABLES := $(SOURCES:.cpp=)

ALL: $(EXECUTABLES)

% : %.cpp $(HEADERS)
	g++ -O3 -std=c++23 -Wall -Wextra -Wshadow -Wpedantic -I.. -o $@ $<

clean:
	rm -f $(EXECUTABLES)

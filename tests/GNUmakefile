SOURCES := $(wildcard *.cpp)
EXECUTABLES := $(SOURCES:.cpp=)

ALL: $(EXECUTABLES)



HEADERS := $(wildcard *.H)

% : %.cpp $(HEADERS)
	g++ -std=c++23 -Wall -Wextra -Wshadow -Wpedantic -I.. -o $@ $<

clean:
	rm -f $(EXECUTABLES)

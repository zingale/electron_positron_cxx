ALL: test_fermi

HEADERS := $(wildcard *.H)
SOURCES := $(wildcard *.cpp)
OBJECTS := $(OBJECTS:.cpp=.o)

%.o : %.cpp $(HEADERS)
	g++ -std=c++23 -c $<

test_fermi : test_fermi.o $(HEADERS)
	g++ -o $@ $<

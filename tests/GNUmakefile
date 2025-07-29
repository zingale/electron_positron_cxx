ALL: test_breakpoints \
     test_brent \
     test_difference \
     test_fermi



HEADERS := $(wildcard *.H)

% : %.cpp $(HEADERS)
	g++ -std=c++23 -Wall -Wextra -Wshadow -I.. -o $@ $<

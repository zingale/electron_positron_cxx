HEADERS := $(wildcard ../*.H)
SOURCES := $(wildcard *.cpp)
EXECUTABLES := $(SOURCES:.cpp=)

#DEFINES := -DUSE_LONG_DOUBLE
DEFINES := -DUSE_FLOAT128

INCLUDES := ..

INCLUDE_FLAGS := $(patsubst %,-I%,$(INCLUDES))
CXXFLAGS := -O3 -std=c++23 -Wall -Wextra -Wshadow

GNU_QUADMATH_FLAGS := -fext-numeric-literals 

# for clang-tidy
PATH_TO_QUADMATH := /usr/lib/gcc/x86_64-redhat-linux/15/include

ALL: $(EXECUTABLES)

% : %.cpp $(HEADERS)
ifeq ($(USE_CLANG_TIDY),TRUE)
	clang-tidy --config-file=../.clang-tidy $< -- $(CXXFLAGS) $(INCLUDE_FLAGS) -I$(PATH_TO_QUADMATH)
endif
	g++  $(CXXFLAGS) $(GNU_QUADMATH_FLAGS) $(DEFINES) $(INCLUDE_FLAGS) -o $@ $<

clean:
	rm -f $(EXECUTABLES)

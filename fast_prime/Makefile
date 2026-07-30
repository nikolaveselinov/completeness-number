CXX ?= g++
CXXFLAGS ?= -O3 -std=c++20 -DNDEBUG -Wall -Wextra -Wpedantic

.PHONY: all clean

all: prime_norm_counter

prime_norm_counter: prime_norm_counter.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	$(RM) prime_norm_counter

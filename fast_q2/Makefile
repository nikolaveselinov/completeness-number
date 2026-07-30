CXX ?= g++
CXXFLAGS ?= -O3 -std=c++20 -DNDEBUG -Wall -Wextra -Wpedantic

.PHONY: all clean

all: q2_norm_counter q2_walsh_counter

q2_norm_counter: q2_norm_counter.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

q2_walsh_counter: q2_walsh_counter.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	$(RM) q2_norm_counter q2_walsh_counter

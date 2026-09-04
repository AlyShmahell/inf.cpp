# CXX=g++ or clang++ — C++17
CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

.PHONY: all compile unit rsa clean

all: unit rsa

compile: unit-bin rsa-bin

unit-bin: unit.cpp ../inf.hpp
	$(CXX) $(CXXFLAGS) -o unit unit.cpp

rsa-bin: rsa.cpp ../inf.hpp
	$(CXX) $(CXXFLAGS) -o rsa rsa.cpp

unit: unit-bin
	./unit

rsa: rsa-bin
	./rsa

clean:
	rm -f unit rsa

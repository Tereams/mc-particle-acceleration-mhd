CXX = g++
CXXFLAGS = -O2 -std=c++17 -Wall -I./src

SRC = main.cpp
TARGET = sim

all:
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET)

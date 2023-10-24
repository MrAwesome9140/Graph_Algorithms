CXX = g++
CXXFLAGS = -std=c++11

SRCS = graph_algorithm.cpp
OUT = graph_algorithm

all: $(OUT)

$(OUT): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $(OUT) $(SRCS)

clean:
	rm -f $(OUT)
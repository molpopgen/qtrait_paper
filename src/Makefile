CXXFLAGS=-O2 -std=c++11

all: summstats.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o summstats summstats.o $(LDFLAGS) -lsequence -lz -lboost_iostreams -ltbb

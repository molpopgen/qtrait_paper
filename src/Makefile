CXXFLAGS=-O2 -std=c++11

all: summstats.o trimFixations.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o summstats summstats.o $(LDFLAGS) -lsequence -lz -lboost_iostreams -ltbb
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o trimFixations trimFixations.o $(LDFLAGS) -lsequence -lz -lboost_iostreams

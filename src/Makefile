CXXFLAGS=-O2 -std=c++11

all: summstats.o trimFixations.o sfs.o tosweed.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o summstats summstats.o $(LDFLAGS) -lsequence -lz -lboost_iostreams -ltbb
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o sfs sfs.o $(LDFLAGS) -lsequence -lz -lboost_iostreams -ltbb
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o trimFixations trimFixations.o $(LDFLAGS) -lsequence -lz -lboost_iostreams
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o tosweed tosweed.o $(LDFLAGS) -lsequence -lz -ltbb

clean: 
	rm -f *.o

#include <Sequence/PolyTableFunctions.hpp>
#include <Sequence/SimData.hpp>
#include <algorithm>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <vector>
using namespace std;
using namespace Sequence;
using namespace boost::iostreams;



int
main(int argc, char **argv)
{
    int argn = 1;
    const char *infile = argv[argn++];

    filtering_istream in;
    ifstream input(infile, ios_base::in | ios_base::binary);
    in.push(gzip_decompressor());
    in.push(input);
    do
        {
            SimData d;
            in >> d >> ws;
	    cout << removeInvariantPos(d) << '\n';
        }
    while (!in.eof());
}

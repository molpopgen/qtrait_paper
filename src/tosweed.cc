/* Convert "big" simulation output from Tennessen model
 * to input files for SweeD
 */

#include <Sequence/SimData.hpp>
#include <Sequence/PolyTableSlice.hpp>
#include <zlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tbb/task_group.h>
#include <vector>
using namespace Sequence;
using namespace std;

void
process_slice(SimData &&d_, const char *stub, const int i, const int window)
{
    string outfile(stub);
    outfile += ".locus" + to_string(i) + ".window" + to_string(window) + ".gz";
    SimData d(std::move(d_));
    ostringstream buffer;
    buffer << d << '\n';
    gzFile out=gzopen(outfile.c_str(), "w");
    gzwrite(out, buffer.str().c_str(), buffer.str().size());
    gzclose(out);
}

void
process_data(SimData &&d_, const char *stub, const int i)
{
    SimData d(std::move(d_));
    double start = static_cast<double>(i) + static_cast<double>(i) * 11.;
    double stop = static_cast<double>(i) + static_cast<double>(i) * 11. + 11.;
    PolyTableSlice<SimData> p(d.sbegin(), d.send(), 1.0, 1.0, start, stop);
    int window = 0;
    for (auto itr = p.cbegin(); itr != p.cend(); ++itr, ++window)
        {
            auto slice = p.get_slice(itr);
			process_slice(std::move(slice),stub,i,window);
        }
}

void
read_data(const char *infile, const char *stub)
{
    std::ifstream in(infile);
    tbb::task_group g;
    SimData d;
    int i = 0;
    while (!in.eof())
        {
            in >> d >> ws;
            g.run([&]() { process_data(std::move(d), stub, i++); });
        }
    g.wait();
}

int
main(int argc, char **argv)
{
    int argn = 1;
    const char *infile = argv[argn++];
    const char *stub = argv[argn++];
    read_data(infile, stub);
}

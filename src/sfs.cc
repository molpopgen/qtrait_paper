/* Calculate "genome-scan" statistics in parallel from
 * files of "ms"-format blocks"
 */
#include <Sequence/PolySIM.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <Sequence/SimData.hpp>
#include <Sequence/SummStats.hpp>
#include <algorithm>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <vector>

using namespace std;
using namespace Sequence;
using namespace boost::iostreams;

// This is our task object,
// which keeps the data,
// and operator() calcuates the statistics
// when desired.
struct sfsTask
{
    string infile;
    explicit sfsTask(string infile_) : infile(infile_) {}
    void
    operator()(void)
    {
        // Read the input file (using boost's iostreams
        // library), and populate a vector of tasks.
        filtering_istream in;
        ifstream input(infile, ios_base::in | ios_base::binary);
        in.push(gzip_decompressor());
        in.push(input);
        // Write the results to a gzip-compressed file
        filtering_ostream out;
        auto outfile = infile;
		auto pos = outfile.find("gz");
		outfile.replace(pos,2,"sfs.gz");
        ofstream output(outfile, ios_base::out | ios_base::binary);
        out.push(gzip_compressor());
        out.push(output);
        SimData d;
        while (!in.eof())
            {
                in >> d >> ws;
                std::vector<unsigned> sfs;
                d = removeInvariantPos(d);
				sfs.resize(d.size()-1);
                for (auto i = d.sbegin(); i != d.send(); ++i)
                    {
                        auto c = std::count(i->second.begin(), i->second.end(),
                                            '1');
                        sfs[c - 1]++;
                    }
                copy(sfs.begin(), sfs.end() - 1,
                     ostream_iterator<unsigned>(out, " "));
                out << sfs.back() << '\n';
            }
        out.pop();
        out.pop();
    }
};

int
main(int argc, char **argv)
{
    int argn = 1;
    const int nthreads = atoi(argv[argn++]);
    tbb::task_scheduler_init init(nthreads);

    vector<sfsTask> tasks;
    for (; argn < argc; ++argn)
        {
            tasks.emplace_back(string(argv[argn]));
        }
    // Execute the tasks in parallel
    tbb::parallel_for(tbb::blocked_range<size_t>(0, tasks.size()),
                      [&tasks](const tbb::blocked_range<size_t> &r) {
                          for (size_t i = r.begin(); i < r.end(); ++i)
                              tasks[i]();
                      });
}

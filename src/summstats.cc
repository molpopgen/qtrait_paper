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

//This is our task object,
//which keeps the data,
//and operator() calcuates the statistics 
//when desired.
struct statTask
{
    SimData d;
    double minfreq, binsize, thetaw, thetapi, tajd, hprime, H1, H12, H2H1, nSL,
        iHS;
    unsigned nd;
    explicit statTask(SimData &&d_, double m, double b)
        : d(move(d_)), minfreq(m), binsize(b),
          thetaw(numeric_limits<double>::quiet_NaN()),
          thetapi(numeric_limits<double>::quiet_NaN()),
          tajd(numeric_limits<double>::quiet_NaN()),
          hprime(numeric_limits<double>::quiet_NaN()),
          H1(numeric_limits<double>::quiet_NaN()),
          H12(numeric_limits<double>::quiet_NaN()),
          H2H1(numeric_limits<double>::quiet_NaN()),
          nSL(numeric_limits<double>::quiet_NaN()),
          iHS(numeric_limits<double>::quiet_NaN()),
          nd(numeric_limits<unsigned>::max())
    {
    }
    void
    operator()(void)
    {
        d = removeInvariantPos(d);
        PolySIM ad(&d);
        thetapi = ad.ThetaPi();
        thetaw = ad.ThetaW();
        tajd = ad.TajimasD();
        hprime = ad.Hprime(true);
        auto G = H1H12(d);
        H1 = G.H1;
        H12 = G.H12;
        H2H1 = G.H2H1;
        nd = ad.NumExternalMutations();
        auto x = snSL(d, minfreq, binsize);
        nSL = x.first;
        iHS = x.second;
    }
    std::ostream &
    write(ostream &o) const
    {
        o << thetaw << ' ' << thetapi << ' ' << nd << ' ' << tajd << ' '
          << hprime << ' ' << H1 << ' ' << H12 << ' ' << H2H1 << ' ' << nSL
          << ' ' << iHS;
        return o;
    }
};

//Our task object is writable to
//and ostream-like type.
std::ostream &
operator<<(ostream &out, const statTask &s)
{
    return s.write(out);
}

int
main(int argc, char **argv)
{
    int argn = 1;
    const char *infile = argv[argn++];
    const char *outfile = argv[argn++];
    const double minfreq = atof(argv[argn++]);
    const double binsize = atof(argv[argn++]);
    const int nthreads = atoi(argv[argn++]);
    tbb::task_scheduler_init init(nthreads);

	//Read the input file (using boost's iostreams
	//library), and populate a vector of tasks.
    filtering_istream in;
    ifstream input(infile, ios_base::in | ios_base::binary);
    in.push(gzip_decompressor());
    in.push(input);
    vector<statTask> tasks;
    do
        {
            SimData d;
            in >> d >> ws;
            tasks.emplace_back(statTask(move(d), minfreq, binsize));
        }
    while (!in.eof());

	//Execute the tasks in parallel
    tbb::parallel_for(tbb::blocked_range<size_t>(0, tasks.size()),
                      [&tasks](const tbb::blocked_range<size_t> &r) {
                          for (size_t i = r.begin(); i < r.end(); ++i)
                              tasks[i]();
                      });
	//Write the results to a gzip-compressed file
    filtering_ostream out;
    ofstream output(outfile, ios_base::out | ios_base::binary);
    out.push(gzip_compressor());
    out.push(output);
    out << "thetaw" << ' ' << "thetapi" << ' ' << "nd" << ' ' << "tajd" << ' '
        << "hprime" << ' ' << "H1" << ' ' << "H12" << ' ' << "H2H1" << ' ' << "nSL"
        << ' ' << "iHS" << '\n';
    copy(tasks.begin(), tasks.end(), ostream_iterator<statTask>(out, "\n"));
    out.pop();
    out.pop();
}

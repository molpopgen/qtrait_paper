#ifndef WRITEMS_HPP
#define WRITEMS_HPP
#include <Sequence/SimData.hpp>
#include <sstream>
#include <zlib.h>
namespace writems {
void writeSimData(const Sequence::SimData & d, const std::string & filename) {
	std::ostringstream o;
	o << d << '\n';
	gzFile gz = gzopen(filename.c_str(),"ab");
	gzwrite(gz,o.str().c_str(),o.str().size());
	gzclose(gz);
}
}

#endif


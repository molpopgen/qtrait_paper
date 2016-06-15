from cython.parallel import parallel, prange
from libsequence.polysitevector cimport polySiteVector,psite_vec_itr,psite_vec_const_itr

from libsequence.summstats cimport PolySIM,GarudStats,H1H12,snSL
from libsequence.polytable cimport SimData,PolyTable
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.utility cimport pair
from libcpp.string cimport string as cppstring
from libcpp.memory cimport unique_ptr
from cython.operator cimport dereference as deref

cdef class simDataVec:
    """
    Wrapper for std::vector<Sequence::SimData>    
    """
    cdef vector[SimData] vec
    def __cinit__(self):
        self.vec = vector[SimData]()
    def __dealloc__(self):
        self.vec.clear()
    def __init__(self,const vector[polySiteVector] & p):
        cdef int i=0
        cdef int n = p.size()
        cdef vector[pair[double,cppstring]].iterator a,b
        while i<n:
            a=p[i].begin()
            b=p[i].end()
            self.vec.push_back(SimData(a,b))
            i+=1
        
cdef map[cppstring,double] get_stats_details(const SimData & d, double minfreq, double binsize ) nogil:
    cdef map[cppstring,double] rv
    cdef unique_ptr[PolySIM] p = unique_ptr[PolySIM](new PolySIM(&d))
    cdef int t = 1
    cdef GarudStats g = H1H12(d)
    cdef double * gmap = NULL
    cdef pair[double,double] nSLiHS = snSL(d,minfreq,binsize,gmap)
    rv[cppstring(b'tajd')] = p.get().TajimasD()
    rv[cppstring(b'thetaw')] = p.get().ThetaW()
    rv[cppstring(b'thetapi')] = p.get().ThetaPi()
    rv[cppstring(b'nd1')] = p.get().NumExternalMutations()
    rv[cppstring(b'hprime')] = p.get().Hprime(t)
    rv[cppstring(b'nSL')] = nSLiHS.first
    rv[cppstring(b'iHS')] = nSLiHS.second
    rv[cppstring(b'H12')] = g.H12
    rv[cppstring(b'H1')] = g.H1
    rv[cppstring(b'H2H1')] = g.H2H1
    
    return rv

cdef vector[map[cppstring,double]] get_stats_parallel( const vector[SimData] & v,double minfreq, double binsize ) nogil:
    cdef vector[map[cppstring,double]] rv
    rv.resize(v.size())
    if v.empty():
        return rv
    cdef int i
    cdef int n = v.size()
    for i in prange(n,schedule='static',nogil=True,chunksize=1):
        rv[i] = get_stats_details(v[i],minfreq,binsize)
    return rv
        
def getSummStatsParallel( simDataVec p, double minfreq = 0.05, double binsize = 0.10 ):
    """
    For the length of p, use that many threads to get summary stats.

    The return value is a list of dicts.  Output order = input order
    """
    if isinstance(simDataVec,p):
        return get_stats_parallel(p.vec,minfreq,binsize)

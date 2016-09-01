from cython_gsl cimport gsl_rng,gsl_rng_get,gsl_rng_alloc,gsl_rng_mt19937,gsl_rng_free
from fwdpy.fwdpp cimport sep_sample_t
from fwdpy.fwdpy cimport TemporalSampler,GSLrng,sampler_base,custom_sampler_data,multilocus_t,uint,sample_sep_single_mloc
from libsequence.polysitevector cimport polySiteVector,psite_vec_itr,psite_vec_const_itr
from libsequence.summstats cimport PolySIM,GarudStats,H1H12,snSL
from libsequence.polytable cimport SimData
from libsequence.extensions cimport SimDataVec
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.utility cimport pair
from libcpp.string cimport string as cppstring
from libcpp.memory cimport unique_ptr
from cython.operator cimport dereference as deref

#This is what gets calculated every generation
ctypedef vector[map[cppstring,double]] final_t

ctypedef pair[uint,gsl_rng *] data_t

#This is the C++ type of our custom sampler
ctypedef custom_sampler_data[final_t,data_t] MlocusSampler_t

#This is the function that actually calculated the summary stats each generation
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

#This is the function that is the workhorse of the sampler
#note: we treat nsam as if const, but cannot declare it const b/c custom sampler API requires the flexibility
cdef void mlocus_sampler_details(const multilocus_t * pop, const unsigned generation, final_t & f, data_t & data) nogil:
    cdef vector[pair[double,double]] b
    cdef vector[sep_sample_t] sample = sample_sep_single_mloc[multilocus_t](data.second,deref(pop),data.first,True,b)
    cdef SimData d
    cdef map[cppstring,double] temp
    for i in range(sample.size()):
        d=SimData(sample[i].first.begin(),sample[i].first.end())
        temp = get_stats_details(d,0.05,0.10)
        temp[cppstring(b'generation')]=<double>generation
        temp[cppstring(b'locus')]=<double>i
        f.push_back(temp)

#Finally, our extension class
cdef class MlocusSummStatsSampler(TemporalSampler):
    def __cinit__(self,unsigned n, unsigned nsam,GSLrng r):
        for i in range(n):
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[MlocusSampler_t](new MlocusSampler_t(&mlocus_sampler_details,
                                                                                                         data_t(nsam,
                                                                                                                gsl_rng_alloc(gsl_rng_mt19937)))))
    def __dealloc__(self):
        """
        Gotta free the gsl_rng *!!!!
        """
        for i in range(self.vec.size()):
            gsl_rng_free((<MlocusSampler_t*>self.vec[i].get()).data.second)

    def get(self):
        cdef vector[final_t] rv
        for i in range(self.vec.size()):
            rv.push_back((<MlocusSampler_t*>self.vec[i].get()).final())
        return rv

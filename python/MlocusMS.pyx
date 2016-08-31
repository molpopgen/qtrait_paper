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

cdef struct sampling_details:
    unsigned nsam;
    cppstring neutral
    cppstring selected

#This is what gets calculated every generation
ctypedef vector[unsigned] final_t

ctypedef pair[sampling_details,gsl_rng *] data_t

#This is the C++ type of our custom sampler
ctypedef custom_sampler_data[final_t,data_t] MlocusMS_t

cdef extern from "writems.hpp" namespace "writems" nogil:
    void writeSimData(const SimData & d, const cppstring & filename) nogil

#This is the function that is the workhorse of the sampler
#note: we treat nsam as if const, but cannot declare it const b/c custom sampler API requires the flexibility
cdef void mlocus_mswriter_details(const multilocus_t * pop, const unsigned generation, final_t & f, data_t & data) nogil:
    f.push_back(generation)
    #Will keep fixations...
    cdef vector[sep_sample_t] sample = sample_sep_single_mloc[multilocus_t](data.second,deref(pop),data.first.nsam,False)
    cdef SimData d
    for i in sample:
        #neutral variants:
        d=SimData(i.first.begin(),i.first.end())
        writeSimData(d,data.first.neutral)
        #non-neutral variants:
        d=SimData(i.second.begin(),i.second.end())
        writeSimData(d,data.first.selected)
    
#Finally, our extension class
cdef class MlocusMSwriter(TemporalSampler):
    def __cinit__(self,unsigned n, unsigned nsam,cppstring nfile,cppstring sfile,GSLrng r):
        cdef sampling_details temp
        temp.nsam=nsam
        for i in range(n):
            nf=nfile+'.'+str(i)+'.gz'
            sf=sfile+'.'+str(i)+'.gz'
            temp.neutral=nf
            temp.selected=sf
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[MlocusMS_t](new MlocusMS_t(&mlocus_mswriter_details,
                                                                                                         data_t(temp,
                                                                                                                gsl_rng_alloc(gsl_rng_mt19937)))))
    def __dealloc__(self):
        """
        Gotta free the gsl_rng *!!!!
        """
        for i in range(self.vec.size()):
            gsl_rng_free((<MlocusMS_t*>self.vec[i].get()).data.second)

    def get(self):
        cdef vector[final_t] rv
        for i in range(self.vec.size()):
            rv.push_back((<MlocusMS_t*>self.vec[i].get()).final())
        return rv

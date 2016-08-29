#Temporal sampler to look at VG per locus over time
from libcpp.memory cimport unique_ptr
from libcpp.utility cimport pair
from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdio cimport printf
from cython_gsl.gsl_statistics cimport *
from fwdpy.fwdpy cimport TemporalSampler,GSLrng,sampler_base,custom_sampler,multilocus_t,uint,sample_sep_single_mloc

cdef struct statdata:
    unsigned generation
    unsigned locus
    string stat
    double value

ctypedef vector[statdata] final_t

ctypedef custom_sampler[final_t] PopstatsLocus_t

cdef void popstats_locus_details(const multilocus_t * pop,
        const unsigned generation,
        final_t & f) nogil:
    #We need to get the 
    #contribution of each locus
    #to trait value 
    #for each individual
    cdef unsigned nloci = pop.diploids[0].size()
    cdef unsigned N = pop.diploids.size()
    cdef vector[double] G
    G.resize(N*nloci,0.0)
    cdef size_t dip,locus,mut
    for dip in range(pop.diploids.size()):
        for locus in range(pop.diploids[dip].size()):
            for mut in range(pop.gametes[pop.diploids[dip][locus].first].smutations.size()):
                G[locus*<size_t>N+dip] += pop.mutations[pop.gametes[pop.diploids[dip][locus].first].smutations[mut]].s
            for mut in range(pop.gametes[pop.diploids[dip][locus].second].smutations.size()):
                G[locus*<size_t>N+dip] += pop.mutations[pop.gametes[pop.diploids[dip][locus].second].smutations[mut]].s
    #Iterate over each locus and get VG
    #for each locus
    cdef double * data = G.data()
    cdef double VG
    cdef statdata temp
    temp.generation=generation
    temp.stat = string("VG")
    for locus in range(nloci):
        VG = gsl_stats_variance(data+locus*<size_t>N,1,N)
        temp.locus=<unsigned> locus
        temp.value=VG
        f.push_back(temp)

cdef class PopstatsLocus(TemporalSampler):
    def __cinit__(self,unsigned n):
        for i in range(n):
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[PopstatsLocus_t](new
                PopstatsLocus_t(&popstats_locus_details)))
    def get(self):
        rv=[]
        for i in range(self.vec.size()):
            rv.append((<PopstatsLocus_t*>self.vec[i].get()).final())
        return rv

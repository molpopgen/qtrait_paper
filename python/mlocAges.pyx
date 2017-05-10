#TemporalSampler to record mean allele age/locus/generation
from libcpp.utility cimport pair
from libcpp.vector cimport vector  
from libcpp.memory cimport unique_ptr 
from fwdpy.fwdpy cimport TemporalSampler,GSLrng,sampler_base,custom_sampler,multilocus_t,uint,sample_separate_mlocus

cdef struct mlocus_allele_age:
    unsigned generation
    unsigned locus
    double age

#generation,locus,mean allele age
ctypedef vector[mlocus_allele_age] final_t

ctypedef custom_sampler[final_t] MlocusAgeSampler_t

cdef void mlocus_age_sampler_details(const multilocus_t * pop,
        const unsigned generation,
        final_t & f) nogil:
    #Keep track of per-locus data 
    #using fast hash tables
    cdef vector[unsigned] locus_counts,locus_data
    cdef unsigned nloci = pop.diploids[0].size()
    locus_counts.resize(nloci,0)
    locus_data.resize(nloci,0)

    #value that we'll keep writing to
    cdef mlocus_allele_age temp
    temp.generation = generation

    cdef unsigned twoN = 2*pop.N
    cdef unsigned locus
    cdef unsigned age
    cdef size_t i
    for i in range(pop.mutations.size()):
        if pop.mutations[i].neutral == 0:
            if pop.mcounts[i]>0 and pop.mcounts[i]<twoN:
                #variant is segregating
                locus = <unsigned>pop.mutations[i].pos
                locus_counts[locus]+=1
                age = (generation - pop.mutations[i].g + 1)
                locus_data[locus] += age
    for i in range(locus_counts.size()):
        if locus_counts[i]:
            temp.locus = i
            temp.age = <double>locus_data[i]/<double>locus_counts[i]
            f.push_back(temp)

cdef class MlocusAgeSampler(TemporalSampler):
    def __cinit__(self,unsigned n):
        for i in range(n):
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[MlocusAgeSampler_t](new
                MlocusAgeSampler_t(&mlocus_age_sampler_details)))

    def get(self):
        rv=[]
        for i in range(self.vec.size()):
            rv.append((<MlocusAgeSampler_t*>self.vec[i].get()).final())
        return rv

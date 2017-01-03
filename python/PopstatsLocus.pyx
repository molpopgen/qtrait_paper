#Temporal sampler to look at VG per locus over time
from libcpp.memory cimport unique_ptr
from libcpp.utility cimport pair
from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdio cimport printf
from libcpp.limits cimport numeric_limits
from fwdpy.gsl cimport *
from cython_gsl.gsl_statistics cimport *
from cython_gsl.gsl_matrix cimport *
from cython_gsl.gsl_vector cimport *
from cython_gsl.gsl_blas cimport *
from cython_gsl.gsl_linalg cimport * 
from cython_gsl.gsl_math cimport * 
from cython_gsl.gsl_sort cimport * 
from fwdpy.fwdpy cimport TemporalSampler,GSLrng,sampler_base,custom_sampler,multilocus_t,uint,sample_sep_single_mloc

from fwdpy.numeric_gsl cimport sum_of_squares 

cdef extern from "<cmath>" namespace "std" nogil:
    bint isnan(float)
    bint isnan(double)
    bint isnan(long double)

cdef extern from "<algorithm>" namespace "std" nogil:
    void reverse[iter](iter,iter)
    int count[iter,value](iter,iter,value)

cdef extern from "<numeric>" namespace "std" nogil:
    oiter partial_sum[iter,oiter](iter beg,iter end,oiter result)

cdef struct statdata:
    unsigned generation
    unsigned locus
    unsigned rank
    double crsq,VG

ctypedef vector[statdata] final_t

ctypedef custom_sampler[final_t] PopstatsLocus_t

cdef void popstats_locus_details(const multilocus_t * pop,
        const unsigned generation,
        final_t & f,
        const int doVG) nogil:
    #We need to get the 
    #contribution of each locus
    #to trait value 
    #for each individual
    cdef unsigned nloci = pop.diploids[0].size()
    cdef unsigned N = pop.diploids.size()
    cdef size_t dip,locus,mut

    cdef gsl_matrix_ptr_t LOCI
    cdef gsl_vector_ptr_t GVALUES
    LOCI.reset(gsl_matrix_alloc(pop.diploids.size(),pop.diploids[0].size()+1))
    gsl_matrix_set_zero(LOCI.get())
    GVALUES.reset(gsl_vector_alloc(pop.diploids.size()))
    gsl_vector_set_zero(GVALUES.get())
    for dip in range(pop.diploids.size()):
        if doVG:
            gsl_vector_set(GVALUES.get(),dip,pop.diploids[dip][0].g)
        else:
            gsl_vector_set(GVALUES.get(),dip,pop.diploids[dip][0].w)
        for locus in range(pop.diploids[dip].size()):
            for mut in range(pop.gametes[pop.diploids[dip][locus].first].smutations.size()):
                s=pop.mutations[pop.gametes[pop.diploids[dip][locus].first].smutations[mut]].s
                gsl_matrix_set(LOCI.get(),dip,locus+1,gsl_matrix_get(LOCI.get(),dip,locus+1) + s)
            for mut in range(pop.gametes[pop.diploids[dip][locus].second].smutations.size()):
                s=pop.mutations[pop.gametes[pop.diploids[dip][locus].second].smutations[mut]].s
                gsl_matrix_set(LOCI.get(),dip,locus+1,gsl_matrix_get(LOCI.get(),dip,locus+1) + s)

    #First things first: 
    #We need to indexes referring to the loci
    #in descending order of VG
    cdef vector[double] VG
    VG.resize(pop.diploids[0].size(),0.0)

    cdef gsl_vector_view column
    for dip in range(1,LOCI.get().size2):
        column = gsl_matrix_column(LOCI.get(),dip)
        VG[dip-1]=gsl_stats_variance(column.vector.data,column.vector.stride,column.vector.size)

    cdef vector[size_t] VGindexes
    VGindexes.resize(VG.size())
    
    gsl_sort_index(VGindexes.data(),VG.data(),1,VG.size())
    reverse(VGindexes.begin(),VGindexes.end())
    
    cdef int invariant = count(VG.begin(),VG.end(),0.0)
    cdef statdata temp
    temp.generation=generation

    if invariant > 0:
        if invariant == <int>VG.size():
            #if there is no VG in the population,
            #create empty entry for this
            #time point and return
            for locus in range(VG.size()):
                temp.locus=locus
                temp.crsq=0.0
                temp.VG=0.0
                temp.rank=locus+1
                f.push_back(temp)
            return 
        #We cannot include these loci in the calculation,
        #and so we'll change size2 (no. columns) to reflect
        #that some columns have no variation
        #LOCI.get().size2 -= invariant

    #Refill the matrix based on sorted order
    cdef size_t dummy=0
    cdef gsl_matrix * m2=gsl_matrix_alloc(LOCI.get().size1,LOCI.get().size2-invariant)
    cdef gsl_vector_view col = gsl_matrix_column(LOCI.get(),0)
    gsl_matrix_set_col(m2,0,&col.vector)
    #for dip in range(pop.diploids.size()):
    dummy=1
    for locus in VGindexes:
        if VG[locus] != 0.0:
            col = gsl_matrix_column(LOCI.get(),locus+1)
            gsl_matrix_set_col(m2,dummy,&col.vector)
            #for mut in range(pop.gametes[pop.diploids[dip][locus].first].smutations.size()):
            #    s=pop.mutations[pop.gametes[pop.diploids[dip][locus].first].smutations[mut]].s
            #    gsl_matrix_set(LOCI.get(),dip,dummy+1,gsl_matrix_get(LOCI.get(),dip,dummy+1) + s)
            #for mut in range(pop.gametes[pop.diploids[dip][locus].second].smutations.size()):
            #    s=pop.mutations[pop.gametes[pop.diploids[dip][locus].second].smutations[mut]].s
            #    gsl_matrix_set(LOCI.get(),dip,dummy+1,gsl_matrix_get(LOCI.get(),dip,dummy+1) + s)
            dummy+=1
    LOCI.reset(m2)
    ssquares = sum_of_squares(GVALUES.get(),LOCI.get()) 
    if isnan(ssquares.first):
        for locus in range(VG.size()):
            temp.locus=locus
            temp.crsq=numeric_limits[double].quiet_NaN()
            temp.VG=numeric_limits[double].quiet_NaN()
            temp.rank=locus+1
            f.push_back(temp)
            return

    ssquares.second.resize(VG.size(),0.0)
    cdef vector[double] csum
    csum.resize(VG.size(),0.0)
    for dip in range(ssquares.second.size()): ssquares.second[dip]/=ssquares.first
    partial_sum(ssquares.second.begin(),ssquares.second.end(),csum.begin())
    dummy=0
    for locus in range(VGindexes.size()):
        temp.locus=VGindexes[locus]
        temp.crsq=csum[locus]
        temp.VG=VG[VGindexes[locus]]
        temp.rank=locus+1
        f.push_back(temp)

cdef void VG_details_wrapper(const multilocus_t * pop,
        const unsigned generation,
        final_t & f) nogil:
    popstats_locus_details(pop,generation,f,1)

cdef void VW_details_wrapper(const multilocus_t * pop,
        const unsigned generation,
        final_t & f) nogil:
    popstats_locus_details(pop,generation,f,0)

cdef class PopstatsLocus(TemporalSampler):
    """
    Contribution of loci to variance in trait value
    """
    def __cinit__(self,unsigned n):
        for i in range(n):
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[PopstatsLocus_t](new
                PopstatsLocus_t(&VG_details_wrapper)))
    def get(self):
        rv=[]
        for i in range(self.vec.size()):
            rv.append((<PopstatsLocus_t*>self.vec[i].get()).final())
        return rv

cdef class PopstatsLocus2(TemporalSampler):
    """
    Contribution of loci to variance in fitness
    """
    def __cinit__(self,unsigned n):
        for i in range(n):
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[PopstatsLocus_t](new
                PopstatsLocus_t(&VW_details_wrapper)))
    def get(self):
        rv=[]
        for i in range(self.vec.size()):
            rv.append((<PopstatsLocus_t*>self.vec[i].get()).final())
        return rv

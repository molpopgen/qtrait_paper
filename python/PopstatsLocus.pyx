#Temporal sampler to look at VG per locus over time
from libcpp.memory cimport unique_ptr
from libcpp.utility cimport pair
from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdio cimport printf
from fwdpy.gsl cimport *
from cython_gsl.gsl_statistics cimport *
from cython_gsl.gsl_matrix cimport *
from cython_gsl.gsl_vector cimport *
from cython_gsl.gsl_blas cimport *
from cython_gsl.gsl_linalg cimport * 
from cython_gsl.gsl_math cimport * 
from cython_gsl.gsl_sort cimport * 
from fwdpy.fwdpy cimport TemporalSampler,GSLrng,sampler_base,custom_sampler,multilocus_t,uint,sample_sep_single_mloc

cdef extern from "<algorithm>" namespace "std" nogil:
    void reverse[iter](iter,iter)
    int count[iter,value](iter,iter,value)

cdef extern from "<numeric>" namespace "std" nogil:
    oiter partial_sum[iter,oiter](iter beg,iter end,oiter result)

cdef struct statdata:
    unsigned generation
    unsigned locus
    unsigned rank
    string stat
    double value

ctypedef vector[statdata] final_t

ctypedef custom_sampler[final_t] PopstatsLocus_t

cdef pair[double,vector[double]] sum_of_squares(gsl_vector * v,
                                                gsl_matrix * m) nogil:
    cdef gsl_vector_ptr_t TAU,SUMS
    cdef gsl_matrix_ptr_t Q,R
    TAU.reset(gsl_vector_alloc(min(m.size1,m.size2)))
    SUMS.reset(gsl_vector_alloc(m.size1))
    Q.reset(gsl_matrix_alloc(m.size1,m.size1))
    R.reset(gsl_matrix_alloc(m.size1,m.size2))

    gsl_linalg_QR_decomp(m,TAU.get())
    gsl_linalg_QR_unpack(m,TAU.get(),Q.get(),R.get())
    gsl_blas_dgemv(CblasTrans,1.0,Q.get(),v,0.0,SUMS.get())

    cdef pair[double,vector[double]] rv

    cdef size_t i
    for i in range(0,m.size2): 
        rv.second.push_back(gsl_pow_2(gsl_vector_get(SUMS.get(),i+1)))
   
    cdef size_t DF = m.size2-1
    RSS=0.0
    for i in range(DF+1,SUMS.get().size):
        RSS += gsl_pow_2(gsl_vector_get(SUMS.get(),i))

    cdef double sqi=0.
    rv.first = RSS
    for sqi in rv.second: rv.first += sqi
    return rv

cdef void popstats_locus_details(const multilocus_t * pop,
        const unsigned generation,
        final_t & f) nogil:
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
        gsl_vector_set(GVALUES.get(),dip,pop.diploids[dip][0].g)
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
    temp.stat = string("cVG")

    if invariant > 0:
        if invariant == <int>VG.size():
            #if there is no VG in the population,
            #create empty entry for this
            #time point and return
            for locus in range(VG.size()):
                temp.locus=locus
                temp.value=0.0
                temp.rank=locus+1
                f.push_back(temp)
            gsl_matrix_free(LOCI.get())
            gsl_vector_free(GVALUES.get())
            return 
        #We cannot include these loci in the calculation,
        #and so we'll change size2 (no. columns) to reflect
        #that some columns have no variation
        LOCI.get().size2 -= invariant

    #Refill the matrix based on sorted order
    cdef size_t dummy=0
    for dip in range(pop.diploids.size()):
        dummy=0
        for locus in VGindexes:
            if VG[locus] != 0.0:
                for mut in range(pop.gametes[pop.diploids[dip][locus].first].smutations.size()):
                    s=pop.mutations[pop.gametes[pop.diploids[dip][locus].first].smutations[mut]].s
                    gsl_matrix_set(LOCI.get(),dip,dummy+1,gsl_matrix_get(LOCI.get(),dip,dummy+1) + s)
                for mut in range(pop.gametes[pop.diploids[dip][locus].second].smutations.size()):
                    s=pop.mutations[pop.gametes[pop.diploids[dip][locus].second].smutations[mut]].s
                    gsl_matrix_set(LOCI.get(),dip,dummy+1,gsl_matrix_get(LOCI.get(),dip,dummy+1) + s)
                dummy+=1

    ssquares = sum_of_squares(GVALUES.get(),LOCI.get()) 

    ssquares.second.resize(VG.size(),0.0)
    cdef vector[double] csum
    csum.resize(VG.size(),0.0)
    for dip in range(ssquares.second.size()): ssquares.second[dip]/=ssquares.first
    partial_sum(ssquares.second.begin(),ssquares.second.end(),csum.begin())
    dummy=0
    for locus in range(VGindexes.size()):
        temp.locus=VGindexes[locus]
        temp.value=csum[locus]
        temp.rank=locus+1
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

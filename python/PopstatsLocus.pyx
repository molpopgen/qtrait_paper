#Temporal sampler to look at VG per locus over time
from libcpp.memory cimport unique_ptr
from libcpp.utility cimport pair
from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdio cimport printf
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
    cdef size_t dip,locus,mut

    #LOCI is a gsl_matrix *
    LOCI = gsl_matrix_alloc(pop.diploids.size(),pop.diploids[0].size()+1)
    gsl_matrix_set_zero(LOCI)
    GVALUES = gsl_vector_alloc(pop.diploids.size())
    gsl_vector_set_zero(GVALUES)
    for dip in range(pop.diploids.size()):
        gsl_vector_set(GVALUES,dip,pop.diploids[dip][0].g)
        for locus in range(pop.diploids[dip].size()):
            for mut in range(pop.gametes[pop.diploids[dip][locus].first].smutations.size()):
                s=pop.mutations[pop.gametes[pop.diploids[dip][locus].first].smutations[mut]].s
                gsl_matrix_set(LOCI,dip,locus+1,gsl_matrix_get(LOCI,dip,locus+1) + s)
            for mut in range(pop.gametes[pop.diploids[dip][locus].second].smutations.size()):
                s=pop.mutations[pop.gametes[pop.diploids[dip][locus].second].smutations[mut]].s
                gsl_matrix_set(LOCI,dip,locus+1,gsl_matrix_get(LOCI,dip,locus+1) + s)

    #First things first: 
    #We need to indexes referring to the loci
    #in descending order of VG
    cdef vector[double] VG
    VG.resize(pop.diploids[0].size(),0.0)

    cdef gsl_vector_view column
    for dip in range(1,LOCI.size2):
        column = gsl_matrix_column(LOCI,dip)
        VG[dip-1]=gsl_stats_variance(column.vector.data,column.vector.stride,column.vector.size)

    cdef vector[size_t] VGindexes
    VGindexes.resize(VG.size())
    
    gsl_sort_index(VGindexes.data(),VG.data(),1,VG.size())
    reverse(VGindexes.begin(),VGindexes.end())
    
    cdef int invariant = count(VG.begin(),VG.end(),0.0)

    if invariant > 0:
        #We cannot include these loci in the calculation
        gsl_matrix_free(LOCI)
        LOCI=gsl_matrix_alloc(pop.diploids.size(),pop.diploids[0].size()-invariant)
        gsl_matrix_set_zero(LOCI)

    #Refill the matrix based on sorted order
    cdef size_t dummy=0
    for dip in range(pop.diploids.size()):
        dummy=0
        for locus in VGindexes:
            if VG[locus] != 0.0:
                for mut in range(pop.gametes[pop.diploids[dip][locus].first].smutations.size()):
                    s=pop.mutations[pop.gametes[pop.diploids[dip][locus].first].smutations[mut]].s
                    #gsl_matrix_set(LOCI,dip,dummy+1,gsl_matrix_get(LOCI,dip,dummy+1) + s)
                for mut in range(pop.gametes[pop.diploids[dip][locus].second].smutations.size()):
                    s=pop.mutations[pop.gametes[pop.diploids[dip][locus].second].smutations[mut]].s
                    #gsl_matrix_set(LOCI,dip,dummy+1,gsl_matrix_get(LOCI,dip,dummy+1) + s)
                dummy+=1

    

    #We now regress individual genetic values onto the genetic values due to 
    #each locus.  We use QR decomposition for the regression.
    TAU=gsl_vector_alloc(min(LOCI.size1,LOCI.size2))
    SUMS=gsl_vector_alloc(LOCI.size1)
    Q=gsl_matrix_alloc(LOCI.size1,LOCI.size1)
    R=gsl_matrix_alloc(LOCI.size1,LOCI.size2)

    gsl_linalg_QR_decomp(LOCI,TAU)
    gsl_linalg_QR_unpack(LOCI,TAU,Q,R)
    gsl_blas_dgemv(CblasTrans,1.0,Q,GVALUES,0.0,SUMS)

    cdef vector[double] squares

    for dip in range(SUMS.size-1): #re-using dip as dummy here!
        squares.push_back(gsl_pow_2(gsl_vector_get(SUMS,dip+1)))
   
    cdef size_t DF = LOCI.size2-1
    cdef float RSS = 0.0
    for dip in range(DF+1,SUMS.size):
        RSS += gsl_pow_2(gsl_vector_get(SUMS,dip))

    cdef double sqi=0.
    cdef double SS = RSS
    for sqi in squares: SS += sqi

    gsl_matrix_free(LOCI)
    gsl_matrix_free(Q)
    gsl_matrix_free(R)
    gsl_vector_free(GVALUES)
    gsl_vector_free(TAU)
    gsl_vector_free(SUMS)

    dummy=0
    cdef statdata temp
    temp.generation=generation
    temp.stat = string("pVG")
    for locus in range(VGindexes.size()):
        temp.locus=VGindexes[locus]
        if VG[VGindexes[locus]] == 0.:
            temp.value=0.0
        else:
            sqi=0.
            dummy=0
            for dip in range(locus):
                if VG[VGindexes[dip]] != 0.:
                    sqi += squares[dummy]
                    dummy+=1
            sqi /= SS
            temp.value=sqi
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

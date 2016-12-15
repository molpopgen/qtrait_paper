#Calculate various forms of genetic load.
#For revisions of Sanjak et al.
#Load is defined as deviation from optimum 
#fitness.  
#For the Gaussian fitness fxn, that optimum
#value is 1.
#We calculate the following loads:
#Total = 1 - fitness
#Fixed = 1 - (fitness calculated only using fixed variants)
#Segregating = 1 - (fitness calculated only using segregating variants)
#This is a custom plugin based on fwdpy.
#What we return is the mean value of each load over time.

from fwdpy.fwdpy cimport TemporalSampler,sampler_base,custom_sampler_data,singlepop_t,uint
from fwdpy.fitness cimport het_mult_update, hom_mult_update_2
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.memory cimport unique_ptr
from libcpp.unordered_map cimport unordered_map
from libc.math cimport sqrt,pow,exp
from cython.operator cimport dereference as deref

cdef struct load_values:
    unsigned generation
    double total
    double fixed
    double seg
    double total_muts
    double total_fixed
    double total_Aa
    double total_aa

#Functions for re-use below
cdef double gaussian_fitness(double P, double opt,double VS) nogil:
    """
    We always assume opt=0.,VS=1
    """
    return exp(-1.*pow(opt-P,2.)/(2.*VS))

cdef pair[double,double] haplotype_sums_seg(const singlepop_t * pop, const size_t diploid) nogil:
    """
    Returns sum of s for each haplotype, only for seg variants
    """
    cdef pair[double,double] rv
    rv.first=0.
    rv.second=0.
    cdef unsigned twoN=2*pop.diploids.size()
    cdef size_t mut = 0
    for mut in range(pop.gametes[pop.diploids[diploid].first].smutations.size()):
        if pop.mcounts[pop.gametes[pop.diploids[diploid].first].smutations[mut]] < twoN:
            rv.first += pop.mutations[pop.gametes[pop.diploids[diploid].first].smutations[mut]].s
    for mut in range(pop.gametes[pop.diploids[diploid].second].smutations.size()):
        if pop.mcounts[pop.gametes[pop.diploids[diploid].second].smutations[mut]] < twoN:
            rv.second += pop.mutations[pop.gametes[pop.diploids[diploid].second].smutations[mut]].s
    return rv

cdef unordered_map[size_t,int] get_Aa_aa(const singlepop_t * pop, const size_t diploid) nogil:
    cdef unordered_map[size_t,int] genotypes
    cdef unsigned twoN=2*pop.diploids.size()
    for mut in range(pop.gametes[pop.diploids[diploid].first].smutations.size()):
        if pop.mcounts[pop.gametes[pop.diploids[diploid].first].smutations[mut]] < twoN:
            genotypes[pop.gametes[pop.diploids[diploid].first].smutations[mut]]+=1
    for mut in range(pop.gametes[pop.diploids[diploid].second].smutations.size()):
        if pop.mcounts[pop.gametes[pop.diploids[diploid].second].smutations[mut]] < twoN:
            genotypes[pop.gametes[pop.diploids[diploid].second].smutations[mut]]+=1
    return genotypes

cdef double multiplicative_pheno_seg(const singlepop_t * pop, const size_t diploid) nogil:
    """
    Returns calculation of phenotype under multiplicative model, only using segregating variants.
    This function is a bummer--we need a lookup table to deal with het/homo status
    """
    cdef double P=1.
    cdef size_t mut = 0
    #Count # of times each mutation occurs in this diploid:
    genotypes=get_Aa_aa(pop,diploid)
    for i in genotypes:
        if i.second == 1:
            het_mult_update(P,pop.mutations[i.first])
        elif i.second == 2:
            hom_mult_update_2(P,pop.mutations[i.first])
    return 1.-P

cdef double sum_fixed_effects(const singlepop_t * pop) nogil:
    """
    Returns sum of s across fixations
    """
    cdef unsigned twoN = 2*pop.diploids.size()
    #Fixed load
    cdef double ssum = 0.
    for i in range(pop.mutations.size()):
        if pop.mcounts[i]==twoN:
            ssum += pop.mutations[i].s
    return ssum
    
cdef double prod_fixed_effects(const singlepop_t * pop) nogil:
    """
    Returns prod of s across fixations.

    Note: this is the total effect for a diploid, accounting 
    for fact that everyone is a homozygote
    """
    cdef unsigned twoN = 2*pop.diploids.size()
    #Fixed load
    cdef double sprod = 1.
    for i in range(pop.mutations.size()):
        if pop.mcounts[i]==twoN:
            sprod *= (1.+2.*pop.mutations[i].s)
    return sprod

cdef load_values make_return_value(const singlepop_t * pop,unsigned generation) nogil:
    """
    Convenience func to return initialized struct
    """
    cdef load_values rv
    rv.generation=generation
    rv.total=0.0
    rv.fixed=0.0
    rv.seg=0.0
    rv.total_muts=0.0
    rv.total_Aa=0.0
    rv.total_aa=0.0
    cdef unsigned nfixed=0
    cdef unsigned twoN = 2*pop.diploids.size()
    for i in range(pop.mcounts.size()):
        if pop.mcounts[i]==twoN:
            nfixed+=1
    rv.total_fixed=2.*<double>nfixed
    return rv 

cdef void update_number_mutations(const singlepop_t * pop, const size_t diploid,load_values * l) nogil:
    genotypes = get_Aa_aa(pop,diploid)
    cdef unsigned Aa=0,aa=0
    for i in genotypes:
        if i.second==1:
            Aa+=1
        elif i.second==2:
            aa+=1
    l.total_Aa+=<double>Aa
    l.total_aa+=<double>aa
    l.total_muts+=(<double>Aa+2.*<double>aa)

#The following functions calculate the mean load for each of the 3 models.
#Notes:
#The fixed load is constant for everyone, so we just calculate it once.

cdef load_values additive_load(const singlepop_t * pop,const unsigned generation) nogil:
    rv = make_return_value(pop,generation)
    #2*sum_fixed_effects b/c everyone is a homozygote
    #for a fixation
    cdef double fixed_w = gaussian_fitness(2.*sum_fixed_effects(pop),0.,1.)
    rv.fixed=1.-fixed_w

    #Seg and total loads
    cdef pair[double,double] hapsums
    cdef size_t i = 0
    cdef size_t mut=0
    for i in range(pop.diploids.size()):
        rv.total += (1.-pop.diploids[i].w)
        hapsums = haplotype_sums_seg(pop,i)
        rv.seg += (1.-gaussian_fitness(hapsums.first+hapsums.second,0.,1.))
        update_number_mutations(pop,i,&rv)
    rv.total /= <double>pop.diploids.size()
    rv.seg /= <double>pop.diploids.size()
    rv.total_Aa/=<double>pop.diploids.size()
    rv.total_aa/=<double>pop.diploids.size()
    rv.total_muts/=<double>pop.diploids.size()
    return rv;

cdef load_values gbr_load(const singlepop_t * pop,const unsigned generation) nogil:
    rv = make_return_value(pop,generation)
    #For this model, P = sqrt(h1*h2), where the h terms are
    #sum over effect sizes on each haplotype.
    #For fixations, h1=h2, and thus P = sqrt(h2^2)=h1.
    rv.fixed = 1.-gaussian_fitness(sum_fixed_effects(pop),0.,1.)

    #Seg and total loads
    cdef pair[double,double] hapsums
    cdef size_t i = 0
    cdef size_t mut=0
    cdef double P
    for i in range(pop.diploids.size()):
        rv.total += (1.-pop.diploids[i].w)
        hapsums = haplotype_sums_seg(pop,i)
        P=sqrt(hapsums.first*hapsums.second)
        rv.seg += (1.-gaussian_fitness(P,0.,1.))
        update_number_mutations(pop,i,&rv)
    rv.total /= <double>pop.diploids.size()
    rv.seg /= <double>pop.diploids.size()
    rv.total_Aa/=<double>pop.diploids.size()
    rv.total_aa/=<double>pop.diploids.size()
    rv.total_muts/=<double>pop.diploids.size()
    return rv;

cdef load_values multiplicative_load(const singlepop_t * pop,const unsigned generation) nogil:
    rv = make_return_value(pop,generation)
    rv.fixed = 1.-(1.- gaussian_fitness(2.*prod_fixed_effects(pop),0.,1.))

    #Seg and total loads
    cdef size_t i = 0
    cdef double p
    for i in range(pop.diploids.size()):
        rv.total += (1.-pop.diploids[i].w)
        p = multiplicative_pheno_seg(pop,i)
        rv.seg += (1.-gaussian_fitness(p,0.,1.))
        update_number_mutations(pop,i,&rv)
    rv.total /= <double>pop.diploids.size()
    rv.seg /= <double>pop.diploids.size()
    rv.total_Aa/=<double>pop.diploids.size()
    rv.total_aa/=<double>pop.diploids.size()
    rv.total_muts/=<double>pop.diploids.size()
    return rv

#Now, we can construct our custom temporal samplers.
#The extension types are implemented in terms of several
#callback functions: the 3 defined above, plus apply_load_calculator.
ctypedef vector[load_values] final_t
ctypedef load_values(*load_calculator_fxn)(const singlepop_t *, const unsigned) nogil
ctypedef custom_sampler_data[final_t,load_calculator_fxn] load_sampler_t

cdef void apply_load_calculator(const singlepop_t * pop,const unsigned generation, final_t & final, load_calculator_fxn & f) nogil:
    final.push_back(f(pop,generation))

cdef get_data(const vector[unique_ptr[sampler_base]] & vec):
    rv=[]
    for i in range(vec.size()):
        temp=(<load_sampler_t*>vec[i].get()).final()
        rv.append(temp)
    return rv

cdef class additiveLoad(TemporalSampler):
    def __cinit__(self,unsigned n):
        for i in range(n):
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[load_sampler_t](new load_sampler_t(&apply_load_calculator,&additive_load)))
    def get(self):
        return get_data(self.vec)

cdef class gbrLoad(TemporalSampler):
    def __cinit__(self,unsigned n):
        for i in range(n):
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[load_sampler_t](new load_sampler_t(&apply_load_calculator,&gbr_load)))
    def get(self):
        return get_data(self.vec)

cdef class multiplicativeLoad(TemporalSampler):
    def __cinit__(self,unsigned n):
        for i in range(n):
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[load_sampler_t](new load_sampler_t(&apply_load_calculator,&multiplicative_load)))
    def get(self):
        return get_data(self.vec)


OUT=tenLocusPopgenMeans.h5  tenLocusQstatMeans.h5 tenLocusSweepSummaries.h5 \
	NumberLociInvolved.csv \
	LociInvolved.h5

all: $(OUT) 

clean: 
	rm -f $(OUT)

tenLocusPopgenMeans.h5 : ../H2_1.0_OPT_0.1_mu_0.00025_sigmu0.25_10regions_popgen.h5   ../H2_1.0_OPT_0.5_mu_0.00025_sigmu0.25_10regions_popgen.h5   ../H2_1.0_OPT_1_mu_0.00025_sigmu0.25_10regions_popgen.h5  ../H2_1.0_OPT_0.1_mu_0.001_sigmu0.25_10regions_popgen.h5     ../H2_1.0_OPT_0.5_mu_0.001_sigmu0.25_10regions_popgen.h5     ../H2_1.0_OPT_1_mu_0.001_sigmu0.25_10regions_popgen.h5  ../H2_1.0_OPT_0.1_mu_0.005_sigmu0.25_10regions_popgen.h5     ../H2_1.0_OPT_0.5_mu_0.005_sigmu0.25_10regions_popgen.h5     ../H2_1.0_OPT_1_mu_0.005_sigmu0.25_10regions_popgen.h5 
	python tenLocusPopgenMeans.py

tenLocusQstatMeans.h5 : ../H2_1.0_OPT_0.1_mu_0.00025_sigmu0.25_10regions_stats.h5   ../H2_1.0_OPT_0.5_mu_0.00025_sigmu0.25_10regions_stats.h5   ../H2_1.0_OPT_1_mu_0.00025_sigmu0.25_10regions_stats.h5  ../H2_1.0_OPT_0.1_mu_0.001_sigmu0.25_10regions_stats.h5     ../H2_1.0_OPT_0.5_mu_0.001_sigmu0.25_10regions_stats.h5     ../H2_1.0_OPT_1_mu_0.001_sigmu0.25_10regions_stats.h5  ../H2_1.0_OPT_0.1_mu_0.005_sigmu0.25_10regions_stats.h5     ../H2_1.0_OPT_0.5_mu_0.005_sigmu0.25_10regions_stats.h5     ../H2_1.0_OPT_1_mu_0.005_sigmu0.25_10regions_stats.h5 
	#python tenLocusQstatMeans.py

tenLocusSweepSummaries.h5 : ../H2_1.0_OPT_0.1_mu_0.00025_sigmu0.25_10regions_popgen_fixations.h5 ../H2_1.0_OPT_0.1_mu_0.001_sigmu0.25_10regions_popgen_fixations.h5 ../H2_1.0_OPT_0.1_mu_0.005_sigmu0.25_10regions_popgen_fixations.h5 ../H2_1.0_OPT_0.5_mu_0.00025_sigmu0.25_10regions_popgen_fixations.h5 ../H2_1.0_OPT_0.5_mu_0.001_sigmu0.25_10regions_popgen_fixations.h5 ../H2_1.0_OPT_0.5_mu_0.005_sigmu0.25_10regions_popgen_fixations.h5 ../H2_1.0_OPT_1_mu_0.00025_sigmu0.25_10regions_popgen_fixations.h5 ../H2_1.0_OPT_1_mu_0.001_sigmu0.25_10regions_popgen_fixations.h5 ../H2_1.0_OPT_1_mu_0.005_sigmu0.25_10regions_popgen_fixations.h5 tenLocusSweepSummaries.py
	python tenLocusSweepSummaries.py

NumberLociInvolved.csv: LociInvolvedSummaries.py
	python LociInvolvedSummaries.py

LociInvolved.h5: LociInvolvedSummaries.py
	python LociInvolvedSummaries.py

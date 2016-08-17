OUT=tenLocusPopgenMeans.h5  tenLocusQstatMeans.h5 tenLocusSweepSummaries.h5 \
	NumberLociInvolved.csv \
	LociInvolved.h5

all: $(OUT) overallfigsdone.txt images_compressed.txt

clean: 
	rm -f $(OUT) *.tif overallfigsdone.txt images_compressed.txt

POPGENFILES=../H2_1.0_OPT_0.1_mu_0.00025_sigmu0.25_10regions_popgen.h5   ../H2_1.0_OPT_0.5_mu_0.00025_sigmu0.25_10regions_popgen.h5   ../H2_1.0_OPT_1_mu_0.00025_sigmu0.25_10regions_popgen.h5  ../H2_1.0_OPT_0.1_mu_0.001_sigmu0.25_10regions_popgen.h5     ../H2_1.0_OPT_0.5_mu_0.001_sigmu0.25_10regions_popgen.h5     ../H2_1.0_OPT_1_mu_0.001_sigmu0.25_10regions_popgen.h5  ../H2_1.0_OPT_0.1_mu_0.005_sigmu0.25_10regions_popgen.h5     ../H2_1.0_OPT_0.5_mu_0.005_sigmu0.25_10regions_popgen.h5     ../H2_1.0_OPT_1_mu_0.005_sigmu0.25_10regions_popgen.h5 

QTRAITFILES=../H2_1.0_OPT_0.1_mu_0.00025_sigmu0.25_10regions_stats.h5   ../H2_1.0_OPT_0.5_mu_0.00025_sigmu0.25_10regions_stats.h5   ../H2_1.0_OPT_1_mu_0.00025_sigmu0.25_10regions_stats.h5  ../H2_1.0_OPT_0.1_mu_0.001_sigmu0.25_10regions_stats.h5     ../H2_1.0_OPT_0.5_mu_0.001_sigmu0.25_10regions_stats.h5     ../H2_1.0_OPT_1_mu_0.001_sigmu0.25_10regions_stats.h5  ../H2_1.0_OPT_0.1_mu_0.005_sigmu0.25_10regions_stats.h5     ../H2_1.0_OPT_0.5_mu_0.005_sigmu0.25_10regions_stats.h5     ../H2_1.0_OPT_1_mu_0.005_sigmu0.25_10regions_stats.h5 

FIXATIONFILES= ../H2_1.0_OPT_0.1_mu_0.00025_sigmu0.25_10regions_popgen_fixations.h5 ../H2_1.0_OPT_0.1_mu_0.001_sigmu0.25_10regions_popgen_fixations.h5 ../H2_1.0_OPT_0.1_mu_0.005_sigmu0.25_10regions_popgen_fixations.h5 ../H2_1.0_OPT_0.5_mu_0.00025_sigmu0.25_10regions_popgen_fixations.h5 ../H2_1.0_OPT_0.5_mu_0.001_sigmu0.25_10regions_popgen_fixations.h5 ../H2_1.0_OPT_0.5_mu_0.005_sigmu0.25_10regions_popgen_fixations.h5 ../H2_1.0_OPT_1_mu_0.00025_sigmu0.25_10regions_popgen_fixations.h5 ../H2_1.0_OPT_1_mu_0.001_sigmu0.25_10regions_popgen_fixations.h5 ../H2_1.0_OPT_1_mu_0.005_sigmu0.25_10regions_popgen_fixations.h5 tenLocusSweepSummaries.py

tenLocusPopgenMeans.h5 : $(POPGENFILES)
	python tenLocusPopgenMeans.py

tenLocusQstatMeans.h5 : $(QTRAITFILES)
	#python tenLocusQstatMeans.py

tenLocusSweepSummaries.h5 : $(FIXATIONFILES)
	python tenLocusSweepSummaries.py

NumberLociInvolved.csv: LociInvolvedSummaries.py
	python LociInvolvedSummaries.py

LociInvolved.h5: LociInvolvedSummaries.py
	python LociInvolvedSummaries.py

overallfigsdone.txt: SummaryStats.ipynb tenLocusPopgenMeans.h5
	jupyter nbconvert --execute SummaryStats.ipynb
	touch overallfigsdone.txt

images_compressed.txt: overallfigsdone.txt
	bash compress_image.sh
	touch images_compressed.txt

OUT=tenLocusPopgenMeans.h5  tenLocusQstatMeans.h5 tenLocusSweepSummaries.h5 \
	NumberLociInvolved.csv \
	LociInvolved.h5 \
    tenLocusPopgenMeansLociWithFixations.h5 \
    tenLocusPopgenMeansLociWithoutFixations.h5 \

all: $(OUT) overallfigsdone.txt locicontributingplots.txt
	jupyter nbconvert --execute SummaryStats.ipynb

clean: 
	rm -f $(OUT) *.tif overallfigsdone.txt locicontributingplots.txt

compressed:
	ls -1 *.tif | grep -v compressed > tiffiles
	parallel bash compress_image.sh :::: tiffiles
#	find . -name '*.tif' | grep -v compressed | xargs rm -f 
	rm -f tiffiles

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

overallfigsdone.txt: plot_summary_stats.py tenLocusPopgenMeans.h5
	python plot_summary_stats.py tenLocusPopgenMeans.h5 _overall
	touch overallfigsdone.txt

tenLocusPopgenMeansLociWithFixations.h5: tenLocusSweepSummaries.h5 $(POPGENFILES) VariationLociContributing.py loci_contributing.py
	python VariationLociContributing.py

tenLocusPopgenMeansLociWithoutFixations.h5: tenLocusSweepSummaries.h5 $(POPGENFILES) VariationLociNotContributing.py loci_contributing.py
	python VariationLociNotContributing.py

locicontributingplots.txt: tenLocusPopgenMeansLociWithFixations.h5 tenLocusPopgenMeansLociWithoutFixations.h5
	python plot_summary_stats.py tenLocusPopgenMeansLociWithFixations.h5 _with_fixations
	python plot_summary_stats.py tenLocusPopgenMeansLociWithoutFixations.h5 _without_fixations
	touch locicontributingplots.txt


# Polygenic adaptation to optimum shifts

## Validation of simulation (additive model)

Under certain conditions, the expectation of VG is 4*mu*VS.  These sims check that the code is giving that.

To run:

```bash
qsub hpc/validate_additive.sh
#Wait a while
python python process_validate_additive_output.py -i popstats.h5 -o popstats.csv
R --no-save --args popstats.csv popstats.pdf < Rscripts/plot_validate_H2_output.R
```

## TODO

* How does $\sigma_\mu$ affect things?
* I should probably do a "track everything" simulation, so that I can track pop-gen stats for interesting replicates.  Alternately, I should require that the tracking API use a diferernt RNG, which is probably way more efficient.  How, though, w/o breaking things?
* How does linkage affect things?
* TFL2013/GBR model

### The question of linkage

The latter should be done via fwdpp's multi-locus API:

* The total mutation rate is $\mu$
* There are $k$ regions, each with mutation rate such that $\sum_i \mu_i = \mu$.
* Each has a within-region recombination rate
* Each has a within-region neutral mutation rate
* The within-region rates should be set to model a 100kb "human" region, b/c that'll give a lot of neutral SNPs to sample
* Vary linkage amongst regions -- 0.1, 0.25, 0.5
* Track everything--samples and mutation frequncies.

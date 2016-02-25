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
# Work flow for making plots, etc.

**Warning:** I am assuming a machine with 512Gb RAM.  Several of the Makefile targets require > 100Gb and a executed assuming many cores are available.  Make -j xxx is not recommended.  Rather, let each job do its own parallelization.

To execute the workflow:

~~~{sh}
make
~~~

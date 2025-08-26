The program ``generate_maxwell_data.cpp`` loops over T and rho and
checks the 3 Maxwell relations and outputs the relative errors to a
file.  The idea is to run it at different precisions and number of
quadrature points and compare the error.

A script to run some combinations can be run as:

```
./create_data.sh
```

This will build it with various options and run.  Note: this can take
a long time, because of the 256-bit case.

A Jupyter notebook to plot the results and data from a previous run is
in `data/`

# eos

This is a sample driver for the equation of state.  It queries for the
density, temperature, and electron fraction and then outputs the
thermodynamic state and the first- and second-derivatives.

You can also pipe the inputs to the command, e.g.,

```
echo 1.e4 1.e7 0.5 | ./eos
```

will set a density of 1.e4 g/cc, temperature of 1.e7 K, and electron
fraction of 0.5.


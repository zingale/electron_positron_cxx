This has a script that loops over T and rho and checks the 3 Maxwell
relations and outputs the relative errors.  The idea is to run it
at different precisions and compare the error.  For example:

```
make clean
make PRECISION=FLOAT128
./generate_maxwell_data > maxwell_128.txt

make clean
make PRECISION=LONG_DOUBLE
./generate_maxwell_data > maxwell_80.txt

make clean
make PRECISION=DOUBLE
./generate_maxwell_data > maxwell_64.txt
```



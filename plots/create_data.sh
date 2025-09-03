#!/bin/bash

# precision

make clean
make PRECISION=DOUBLE QUAD_PTS=200
./generate_maxwell_data

make clean
make PRECISION=LONG_DOUBLE QUAD_PTS=200
./generate_maxwell_data

make clean
make PRECISION=FLOAT128 QUAD_PTS=200
./generate_maxwell_data

make clean
make PRECISION=BOOST256 QUAD_PTS=200
./generate_maxwell_data

# quadrature points

make clean
make PRECISION=FLOAT128 QUAD_PTS=20
./generate_maxwell_data

make clean
make PRECISION=FLOAT128 QUAD_PTS=50
./generate_maxwell_data

make clean
make PRECISION=FLOAT128 QUAD_PTS=100
./generate_maxwell_data

make clean
make PRECISION=FLOAT128 QUAD_PTS=200
./generate_maxwell_data

make clean
make PRECISION=FLOAT128 QUAD_PTS=400
./generate_maxwell_data

make clean
make PRECISION=FLOAT128 QUAD_PTS=800
./generate_maxwell_data

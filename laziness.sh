#!/bin/bash

for((i = 0; i < 8; i+=1));
do
    ./a.out < "testing/maximization/input/input${i}" > "testing/maximization/output/output${i}"
done;

for((i = 0; i < 4; i+=1));
do
    ./a.out < "testing/minimization/input/input${i}" > "testing/minimization/output/output${i}"
done;
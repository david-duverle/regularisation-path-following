#!/bin/bash

cd lcm53_patched/
make
mv lcm ../
cd ../glmpath/
rm *.so
rm *.o
R CMD SHLIB asa_cg.c solve_coxpath.c -o ../solve_coxpath.so
rm *.o
cd ../itemset-ratio/
rm *.so
rm *.o
R CMD SHLIB explore_itemset.c -o ../explore_itemset.so
rm *.o
cd ..

echo "***** BUILD COMPLETE ****"
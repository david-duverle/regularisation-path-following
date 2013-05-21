rm *.so
rm *.o
R CMD SHLIB asa_cg.c solve_coxpath.c -o solve_coxpath.so

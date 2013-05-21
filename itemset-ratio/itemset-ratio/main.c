//
//  main.c
//  itemset-ratio
//
//  Created by Dave on 11/13/12.
//  Copyright (c) 2012 Dave. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mach/vm_statistics.h>
#include <mach/mach_types.h>
#include <mach/mach_init.h>
#include <mach/mach_host.h>


char *load_bool_file(const char *filename, int nrow, int ncol);
double *load_float_file(const char *filename, int nrow, int ncol);
int *load_int_file(const char *filename, int nrow, int ncol);
void get_used_mem();

void get_min_itemset(int *_tot_depth, int *exprs_col_major, int *_nrow, int *_ncol, double *w1, double *w2, int *active, int *_active_len, double *_lambda, double *curmin, int *min_idxs, int *trace);


int main(int argc, const char * argv[])
{
    get_used_mem();
    
    printf("Loading input files\n");
    
    int depth = 1;
    int active_len = 4;
    double lambda = 3.73749564429404;

    int nrow = 137, ncol = 7;
    
    int min_idxs[depth];
    double res_min = lambda;

    char *exprs_row_major = load_bool_file("/Users/dave/Academia/CBRC/SyntheticLethality/Cox Regression/cox-regression-path/itemset-ratio/data/exprs.txt", nrow, ncol);
    double *w1 = load_float_file("/Users/dave/Academia/CBRC/SyntheticLethality/Cox Regression/cox-regression-path/itemset-ratio/data/w1.txt", nrow, 1);
    double *w2 = load_float_file("/Users/dave/Academia/CBRC/SyntheticLethality/Cox Regression/cox-regression-path/itemset-ratio/data/w2.txt", nrow, 1);
    
    int *exprs_col_major = malloc(ncol * nrow * sizeof(int)); //column-major
    for (int j = 0; j < nrow; j++)
        for (int i = 0; i < ncol; i++)
            exprs_col_major[i*nrow+j] = exprs_row_major[j*ncol+i];

    int *active = load_int_file("/Users/dave/Academia/CBRC/SyntheticLethality/Cox Regression/cox-regression-path/itemset-ratio/data/active.txt", active_len, depth);
//    int *active = NULL;
    
    clock_t begin=clock();
    
    get_used_mem();

    int trace = 3;
    get_min_itemset(&depth, exprs_col_major, &nrow, &ncol, w1, w2, active, &active_len, &lambda, &res_min, (int *) &min_idxs, &trace);
    
    get_used_mem();

    clock_t end=clock();
    double diff = (end - begin)*1000/CLOCKS_PER_SEC;
    printf("********\nExecution took: %.3f s.\nFound min: %lf ", diff/1000, res_min);
    for(int i = 0; i < depth; i++)
        printf("[%d]", min_idxs[i]);
    printf("\n");

    free(exprs_col_major);
    free(w1);
    free(w2);
    free(exprs_row_major);
    free(active);
    
    get_used_mem();
    
    return 0;
}


void get_used_mem() {
    vm_size_t page_size;
    mach_port_t mach_port;
    mach_msg_type_number_t count;
    vm_statistics_data_t vm_stats;
    
    mach_port = mach_host_self();
    count = sizeof(vm_stats) / sizeof(natural_t);
    if (KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
        KERN_SUCCESS == host_statistics(mach_port, HOST_VM_INFO,
                                        (host_info_t)&vm_stats, &count))
    {
//        int64_t free_mem = (int64_t)vm_stats.free_count * (int64_t)page_size;
        
        int64_t used_mem = ((int64_t)vm_stats.active_count +
                       (int64_t)vm_stats.inactive_count +
                       (int64_t)vm_stats.wire_count) *  (int64_t)page_size;
    
        printf("\n#######\nUsed mem: %lldk\n", used_mem/1024);
    }
}

char *load_bool_file(const char *filename, int nrow, int ncol) {
    int val;
    char *array = malloc(nrow*ncol);

    FILE* f = fopen(filename, "r");
    if(!f) {
     printf("Wrong input file: %s", filename);
     exit(-10);
    }

    for (int k = 0; k < nrow*ncol; k++) {
     if (fscanf(f, "%d", &val) > 0)
         array[k] = val;
     else {
         printf("Missing values in file: %s. Got %d expected %d", filename, k-1, nrow*ncol);
         exit(-11);
     }
    }

    fclose(f);

    return array;
}
 
 double *load_float_file(const char *filename, int nrow, int ncol) {
     double *array = malloc(nrow*ncol*sizeof(double));
     double val;
     
     FILE* f = fopen(filename, "r");
     if(!f) {
         printf("Wrong input file: %s", filename);
         exit(-10);
     }
 
    for (int k = 0; k < nrow*ncol; k++) {
        if (fscanf(f, "%lf", &val) > 0)
            array[k] = val;
        else {
            printf("Missing values in file: %s. Got %d expected %d", filename, k-1, nrow*ncol);
            exit(-11);
        }
    }
 
     fclose(f);

     return array;
}

int *load_int_file(const char *filename, int nrow, int ncol) {
    int *array = malloc(nrow*ncol*sizeof(int));
    int val;
    
    FILE* f = fopen(filename, "r");
    if(!f) {
        printf("Wrong input file: %s", filename);
        exit(-10);
    }
    
    for (int k = 0; k < nrow*ncol; k++) {
        if (fscanf(f, "%d", &val) > 0)
            array[k] = val;
        else {
            printf("Missing values in file: %s. Got %d expected %d", filename, k-1, nrow*ncol);
            exit(-11);
        }
    }
    
    fclose(f);
    
    return array;
}

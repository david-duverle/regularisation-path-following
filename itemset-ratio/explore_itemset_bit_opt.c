//
//  main.cpp
//  cox itemset
//
//  Created by Dave on 11/9/12.
//  Copyright (c) 2012 CBRC. All rights reserved.
//

#include <R.h>

// #include <fstream>
// #include <iostream>
// #include <sstream>
// #include <string>
#include <float.h>
#include <stdbool.h>
#include <stdlib.h> 
#include <time.h>
#include <stdint.h>

#include <pthread.h>

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

#define EPSILON 0.0001f

#define MY_MALLOC(x) malloc(x)
//#define MY_MALLOC(x) (void *) R_alloc(x, 1)
//#define MY_MALLOC(x) R_chk_calloc(x, 1)

#define MY_FREE(x) free(x)
//#define MY_FREE(x) // x
//#define MY_FREE(x) Free(x)
#define CHECK_PTR(var, msg) if(var == NULL) { printf(msg); exit(-1); }

#define REGISTER_SIZE 64

static inline void set_bit(char *x, long bitNum) {
    x[bitNum/8] |= (1L << (bitNum % 8));
}

double tot_explored;
int g_trace;

struct get_min_arg_struct {
    //read-only:
    int tot_depth;
    char *exprs;
    char *single_exprs;
    int nrow;
    int ncol;
    double *w1p, *w1m, *w2p, *w2m;
    double num_const, denom_const;
    int *active; int active_len;
    //mutualised:
    double *curmin; int *combi_min_idx;
    pthread_mutex_t *mutex;
};

void get_min_recur(int cur_depth, int tot_depth, int *cur_path, 
                   char *cur_exprs, char *single_exprs, int nrow, int ncol,
                   double *w1p, double *w1m, double *w2p, double *w2m,
                   double num_const, double denom_const,
                   double *curmin, int *combi_min_idx, int *active, int active_len, pthread_mutex_t *mutex, double *mem_buffer);
static inline bool is_idx_in_set(int *cur_path, int tot_depth, int min_idx, int *active, int active_len);

void launch_get_min_recur(void *args);

/*
double get_max(double *arr, int len) {
    double maxval = -DBL_MAX;
    for (int i = 0; i < len; i++)
        maxval = max(arr[i], maxval);
    
    return maxval;
}
*/

static inline int intcmp(const void *aa, const void *bb)
{
    const int *a = aa, *b = bb;
    return (*a < *b) ? -1 : (*a > *b);
}


 int get_pos_min_idx(double *arr, int len, int *cur_path, int tot_depth, int *active, int active_len) {
    double minval = DBL_MAX;
    int min_idx = -1;
    
    for (int i = 0; i < len; i++) {
        if(arr[i] > EPSILON && arr[i] < minval && !is_idx_in_set(cur_path, tot_depth, i, active, active_len)) {
            minval = arr[i];
            min_idx = i;
        }
    }
    return min_idx;
}


void get_min_itemset(int *_tot_depth, int *exprs_col_major, int *_nrow, int *_ncol, double *w1, double *w2, int *active, int *_active_len, double *_lambda, double *curmin, int *min_idxs, int *trace) {
    
    g_trace = *trace;
    
   int tot_depth = *_tot_depth,
      nrow = *_nrow,
      ncol = *_ncol,
      active_len = *_active_len;
   double lambda = *_lambda;
   
    double *w1p = MY_MALLOC(nrow * sizeof(double));
    CHECK_PTR(w1p, "w1p could not be allocated");
    double *w1m = MY_MALLOC(nrow * sizeof(double));
    CHECK_PTR(w1m, "w1m could not be allocated");
    double *w2p = MY_MALLOC(nrow * sizeof(double));
    CHECK_PTR(w2p, "w2p could not be allocated");
    double *w2m = MY_MALLOC(nrow * sizeof(double));
    CHECK_PTR(w2m, "w2m could not be allocated");
    
    int nrow_padded = (nrow / REGISTER_SIZE + 1) * REGISTER_SIZE; // padded to 64 bits
    char *exprs = MY_MALLOC(ncol * nrow_padded / CHAR_BIT); //column-major
    CHECK_PTR(exprs, "exprs could not be allocated");
    memset(exprs, ncol * nrow_padded / CHAR_BIT, 0);
    
    for (int i = 0; i < nrow; i++) {
        w1p[i] = max(0., w1[i]);
        w1m[i] = max(0., -w1[i]);
        w2p[i] = max(0., w2[i]);
        w2m[i] = max(0., -w2[i]);
    }
    
    for (int i = 0; i < active_len*tot_depth; i++)
        active[i]--;
    
    for (int i = 0; i < tot_depth; i++)
        min_idxs[i] = -1;
    

    for (int j = 0; j < nrow; j++)
        for (int i = 0; i < ncol; i++)
            if(exprs_col_major[i*nrow+j])
                set_bit(exprs, i*nrow_padded+j);
            // exprs[i*nrow+j] = exprs_row_major[j*ncol+i];
            
    if(*curmin == 0 || *curmin > lambda)
        *curmin = lambda;
        
    tot_explored = 0;
    
    
    struct get_min_arg_struct args;
    args.tot_depth = tot_depth;
    args.exprs = exprs;
    args.single_exprs = MY_MALLOC(ncol * nrow_padded / CHAR_BIT);
    CHECK_PTR(args.single_exprs, "args.single_exprs could not be allocated");

    memcpy(args.single_exprs, exprs, ncol * nrow_padded / CHAR_BIT);
    args.nrow = nrow;
    args.ncol = ncol;
    args.w1p = w1p;
    args.w1m = w1m;
    args.w2p = w2p;
    args.w2m = w2m;
    args.num_const = lambda;
    args.denom_const = 1;
    args.curmin = curmin;
    args.combi_min_idx = min_idxs;
    args.active = active;
    args.active_len = active_len;
    
    pthread_mutex_t mutex;
    args.mutex = &mutex;
    
    pthread_mutex_init(args.mutex, NULL);

    pthread_t sub1, sub2;
    
    pthread_create(&sub1, NULL, (void *) launch_get_min_recur, &args);
    
    struct get_min_arg_struct args2 = args;
    args2.num_const = -lambda;
    args2.denom_const = -1;
    args2.single_exprs = MY_MALLOC(ncol * nrow_padded / CHAR_BIT);
    CHECK_PTR(args2.single_exprs, "args2.single_exprs could not be allocated");

    memcpy(args2.single_exprs, exprs, ncol * nrow_padded / CHAR_BIT);

    pthread_create(&sub2, NULL, (void *) launch_get_min_recur, &args2);
    
    pthread_join(sub1, NULL);
    pthread_join(sub2, NULL);
    
    pthread_mutex_destroy(args.mutex);



    for (int i = 0; i < args.tot_depth; i++) {
      min_idxs[i] = min_idxs[i]+1;  //IMPORTANT: must increment!
   }
    
    if(g_trace >= 2)
        printf("Total explored: %.0f\n\n", tot_explored);

    MY_FREE(w1p);
    MY_FREE(w1m);
    MY_FREE(w2p);
    MY_FREE(w2m);
    MY_FREE(exprs);
   
   // printf(" (lambda: %.4f)\n", lambda);
}

void launch_get_min_recur(void *args) {
    struct get_min_arg_struct *arg_struct = (struct get_min_arg_struct *)args;
    int *cur_path = MY_MALLOC(arg_struct->tot_depth * sizeof(int));
    CHECK_PTR(cur_path, "cur_path could not be allocated");
    
    for (int i = 0; i < arg_struct->tot_depth; i++)
        cur_path[i] = -1;
    
    double *mem_buffer = MY_MALLOC(5 * arg_struct->ncol * sizeof(double));
    CHECK_PTR(mem_buffer, "mem_buffer could not be allocated");

    if(g_trace >= 2)
        printf("Launching search for (%lf,%lf)\n", arg_struct->num_const, arg_struct->denom_const);
    
    get_min_recur(0, arg_struct->tot_depth, cur_path, arg_struct->exprs, arg_struct->single_exprs, arg_struct->nrow, arg_struct->ncol, arg_struct->w1p, arg_struct->w1m, arg_struct->w2p, arg_struct->w2m, arg_struct->num_const, arg_struct->denom_const, arg_struct->curmin, arg_struct->combi_min_idx, arg_struct->active, arg_struct->active_len, arg_struct->mutex, mem_buffer);
    
    MY_FREE(mem_buffer);
    if(g_trace >= 2)
        printf("Finished search for (%lf,%lf)\n", arg_struct->num_const, arg_struct->denom_const);
}

bool is_idx_in_set(int *cur_path, int tot_depth, int min_idx, int *active, int active_len) {
    int *idxs = MY_MALLOC(tot_depth * sizeof(int));
    CHECK_PTR(idxs, "idxs could not be allocated");

    int i, j;
    
    
    for (i = 0; i < tot_depth && cur_path[i] >= 0; i++)
        idxs[i] = cur_path[i];
    idxs[i] = min_idx;
    for (j = i+1; j < tot_depth; j++)
        idxs[j] = -1;
    
    qsort(idxs, i+1, sizeof(int), intcmp);

    for(i = 0; i < active_len; i++) {
        for (j = 0; j < tot_depth; j++) {
            if(active[i*tot_depth + j] != idxs[j])
                break;
        }
        if(j >= tot_depth) {
            MY_FREE(idxs);
            return true;
        }
    }
    
    MY_FREE(idxs);
    return false;
}

void get_min_recur(int cur_depth, int tot_depth, int *cur_path,
                  char *cur_exprs, char *single_exprs, int nrow, int ncol,
                  double *w1p, double *w1m, double *w2p, double *w2m,
                  double num_const, double denom_const,
                  double *curmin, int *combi_min_idx,
                     int *active, int active_len, pthread_mutex_t *mutex, double *mem_buffer) {
    
    tot_explored++;

    
    double *a = mem_buffer;
    double *b = mem_buffer+ncol;
    double *c = mem_buffer + 2*ncol;
    double *d = mem_buffer + 3*ncol;
    double *gamma = mem_buffer + 4*ncol;
    long nth_bit;
    
    int nrow_padded = (nrow / REGISTER_SIZE + 1) * REGISTER_SIZE; // padded to 64 bits

    for (int j = 0; j < ncol; j++) {
        a[j] = 0;
        b[j] = 0;
        c[j] = 0;
        d[j] = 0;
        for (int i = 0; i < nrow; i++) {
            nth_bit = i + j * nrow_padded;
            if(cur_exprs[nth_bit/8] & (1 << nth_bit%8)) {
//            if(cur_exprs[i + j*nrow]) {
                a[j] += w1p[i];
                b[j] += w1m[i];
                c[j] += w2p[i];
                d[j] += w2m[i];
            }
        }
    }
        
    for (int j = 0; j < ncol; j++)
        gamma[j] = (num_const + a[j] - b[j]) / (denom_const + c[j] - d[j]);
         
    int min_idx = get_pos_min_idx(gamma, ncol, cur_path, tot_depth, active, active_len);

    if (min_idx >= 0 && gamma[min_idx] < *curmin) {
        bool is_duplicate = false;
        for(int i = 0; i < cur_depth; i++)
            if(cur_path[i] == min_idx) {
                is_duplicate = true;
                break;
            }
        
        pthread_mutex_lock(mutex);
        if(! is_duplicate && gamma[min_idx] < *curmin) { // Checking again here to avoid race conditions
            
            *curmin = gamma[min_idx];
            for (int i = 0; i < tot_depth; i++)
                combi_min_idx[i] = cur_path[i];
            combi_min_idx[cur_depth] = min_idx;
            
            if(g_trace >= 3) {
                printf("New min: %.3f ", *curmin);
                for (int i = 0; i < tot_depth; i++)
                   printf("[%d]", combi_min_idx[i]);
                printf("\n");
            }
            
        }
        pthread_mutex_unlock(mutex);
    }
    

    if(cur_depth+1 < tot_depth) {
        
        bool to_skip[ncol];
        for(int j = 0; j < ncol; j++) {
            to_skip[j] = ((denom_const - d[j] >= 0
                 && (((num_const - b[j]) / (denom_const + c[j]) >= *curmin) ||
                     ((num_const + a[j]) / (denom_const - d[j]) <= EPSILON)))
                || (denom_const + c[j] <= 0
                    && (((num_const + a[j]) / (denom_const - d[j]) >= *curmin) ||
                        ((num_const - b[j]) / (denom_const + c[j]) <= EPSILON))));
        }
        
        char *combi_exprs = MY_MALLOC((ncol * nrow / 8) + 1);
        CHECK_PTR(combi_exprs, "combi_exprs could not be allocated");
//        char *exprs_mask = MY_MALLOC((ncol * nrow / 8) + 1);
//        if(exprs_mask == NULL) {
//            printf("exprs_mask could not be allocated: %d * %d / 8 + 1 = %d\ncurdepth: %d\n", ncol, nrow, nrow * (ncol / 8), cur_depth);
//            exit(-1);
//        }
        
        int explored = 0;
        
        clock_t begin=clock();

        char null_block [nrow];
        bzero(null_block, nrow);

        for(int j = 0; j < ncol; j++) {
            for(int i = 0; i < cur_depth; i++)
                if(cur_path[i] == j)
                    goto skip;
            

            if(to_skip[j])
                continue;
            
            explored++;

            bool non_null = false;
            int nrow64 = nrow_padded / REGISTER_SIZE;
            for(int k = 0; k < ncol; k++) {
                for (int i = 0; i < nrow64; i++) {
                    ((uint64_t *) combi_exprs)[k*nrow64 + i] = ((uint64_t *) cur_exprs)[k*nrow64 + i] & ((uint64_t *) single_exprs)[j*nrow64 + i];
                    non_null |= (((uint64_t *) combi_exprs)[k*nrow64 + i] != 0);
                }
            }


            if (! non_null)
                continue;
            
//            for (int i = 0; i < nrow*ncol; i++)
//                if((debug[i] & exprs_mask[i]) != combi_exprs[i])
//                    cout << "Not matching! " << i << endl;
            
            cur_path[cur_depth] = j;
            get_min_recur(cur_depth+1, tot_depth, cur_path, combi_exprs, single_exprs, nrow, ncol, w1p, w1m, w2p, w2m, num_const, denom_const, curmin, combi_min_idx, active, active_len, mutex, mem_buffer);
            
            
            skip:;
            
        }
        cur_path[cur_depth] = -1;
        
        MY_FREE(combi_exprs);
//        MY_FREE(exprs_mask);
        
        clock_t end=clock();
        double diff = (end - begin)*1000/CLOCKS_PER_SEC;
        if(explored > 0 && cur_depth < 1) {
            if(g_trace >= 2)
                printf("Explored %d in %.f ms per branch (total: %.2f s.)\n", explored, diff/explored, diff/1000);
        }
    }
   
    return;
}

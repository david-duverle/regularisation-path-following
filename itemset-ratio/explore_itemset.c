//
//  main.cpp
//  cox itemset
//
//  Created by Dave on 11/9/12.
//  Copyright (c) 2012 CBRC. All rights reserved.
//

#include <R.h>

#include <math.h>
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

#define EPSILON 0.00001f

#ifndef MAXFLOAT
#define MAXFLOAT 10e100
#endif

#define USE_HEURISTIC_1
#define USE_HEURISTIC_2

#define MY_MALLOC(x) malloc(x)
//#define MY_MALLOC(x) (void *) R_alloc(x, 1)
//#define MY_MALLOC(x) R_chk_calloc(x, 1)

#define MY_FREE(x) free(x)
//#define MY_FREE(x) // x
//#define MY_FREE(x) Free(x)
#define CHECK_PTR(var, msg) if(var == NULL) { printf(msg); exit(-1); }


double tot_node_explored, tot_non_leaf_explored;
int g_trace;
static long pruned = 0, not_pruned = 0, easy_pruned = 0;
static double *weight_ratios;

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
                   double *curmin, int *combi_min_idx, bool *skip_single, int *active, int active_len, pthread_mutex_t *mutex);
static bool is_idx_in_set(int *cur_path, int tot_depth, int min_idx, int *active, int active_len);

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

static inline int order_of_double_cmp(const void *aa, const void *bb)
{
    const int *a = aa, *b = bb;
    
    const double da = weight_ratios[*a], db = weight_ratios[*b];
    return (da < db) ? -1 : (da > db);
//    return (db < da) ? -1 : (db > da);
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
    
    // printf("testing input:\ntot_depth: %d\nnrow: %d\nncol: %d\nw1: %.10f\nw2: %.10f\nactive: ", tot_depth, nrow, ncol, w1[0], w2[0]);
    //     printf("##DEBUG: active set: ");
    //     for(int i=(active_len-1) * tot_depth; i < active_len * tot_depth; i++)
    //        printf("%d ", active[i]);
    //    printf("\n\n");
    // printf("\n");
    // printf("\nlambda: %.4f\ncurmin: %.4f\nexprs:\n", lambda, *curmin);
    // for(int i=0; i < 100; i++)
    //    printf("%d ", (int) exprs_row_major[i]);
    // printf("\n");
    
    double *w1p = MY_MALLOC(nrow * sizeof(double));
    CHECK_PTR(w1p, "w1p could not be allocated");
    double *w1m = MY_MALLOC(nrow * sizeof(double));
    CHECK_PTR(w1m, "w1m could not be allocated");
    double *w2p = MY_MALLOC(nrow * sizeof(double));
    CHECK_PTR(w2p, "w2p could not be allocated");
    double *w2m = MY_MALLOC(nrow * sizeof(double));
    CHECK_PTR(w2m, "w2m could not be allocated");
    char *exprs = MY_MALLOC(ncol * nrow); //column-major
    CHECK_PTR(exprs, "exprs could not be allocated");
    weight_ratios = MY_MALLOC(nrow * sizeof(double));
    
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
            exprs[i*nrow+j] = exprs_col_major[i*nrow+j];
    // exprs[i*nrow+j] = exprs_row_major[j*ncol+i];
    
    if(*curmin == 0 || *curmin > lambda)
        *curmin = lambda;
    
    tot_node_explored = 0;
    tot_non_leaf_explored = 0;
    
    
    struct get_min_arg_struct args;
    args.tot_depth = tot_depth;
    args.exprs = exprs;
//    args.single_exprs = MY_MALLOC(ncol * nrow);
//    CHECK_PTR(args.single_exprs, "args.single_exprs could not be allocated");
//    
//    memcpy(args.single_exprs, exprs, ncol*nrow);
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
    
    pthread_create(&sub2, NULL, (void *) launch_get_min_recur, &args2);
    
    pthread_join(sub1, NULL);
    pthread_join(sub2, NULL);
    
    pthread_mutex_destroy(args.mutex);
    
    for (int i = 0; i < args.tot_depth; i++) {
        min_idxs[i] = min_idxs[i]+1;  //IMPORTANT: must increment!
        // printf("[%d]", min_idxs[i]);
    }
    
    if(g_trace >= 2)
        printf("Total explored: %.0f non-terminal (%.0f total)\n Pruned: %ld + %ld(/%ld)\n\n", tot_non_leaf_explored, tot_node_explored, easy_pruned, pruned, pruned+not_pruned);
    
    
//    MY_FREE(args.single_exprs);
//    MY_FREE(args2.single_exprs);
    
    MY_FREE(w1p);
    MY_FREE(w1m);
    MY_FREE(w2p);
    MY_FREE(w2m);
    MY_FREE(exprs);
    MY_FREE(weight_ratios);
    // printf(" (lambda: %.4f)\n", lambda);
    
}

void launch_get_min_recur(void *args) {
    struct get_min_arg_struct *arg_struct = (struct get_min_arg_struct *)args;
    int *cur_path = MY_MALLOC(arg_struct->tot_depth * sizeof(int));
    CHECK_PTR(cur_path, "cur_path could not be allocated");
    
    for (int i = 0; i < arg_struct->tot_depth; i++)
        cur_path[i] = -1;
    
    arg_struct->single_exprs = MY_MALLOC(arg_struct->ncol * arg_struct->nrow);
    CHECK_PTR(arg_struct->single_exprs, "args.single_exprs could not be allocated");
    
    memcpy(arg_struct->single_exprs, arg_struct->exprs, arg_struct->ncol*arg_struct->nrow);

    bool *skip_single = MY_MALLOC(arg_struct->ncol * sizeof(bool));
    memset(skip_single, 0, arg_struct->ncol * sizeof(bool));
    
    if(g_trace >= 2)
        printf("Launching search for (%lf,%lf)\n", arg_struct->num_const, arg_struct->denom_const);
    
    get_min_recur(0, arg_struct->tot_depth, cur_path, arg_struct->exprs, arg_struct->single_exprs, arg_struct->nrow, arg_struct->ncol, arg_struct->w1p, arg_struct->w1m, arg_struct->w2p, arg_struct->w2m, arg_struct->num_const, arg_struct->denom_const, arg_struct->curmin, arg_struct->combi_min_idx, skip_single, arg_struct->active, arg_struct->active_len, arg_struct->mutex);
    
    if(g_trace >= 2)
        printf("Finished search for (%lf,%lf)\n", arg_struct->num_const, arg_struct->denom_const);
    
    MY_FREE(arg_struct->single_exprs);
    MY_FREE(skip_single);
    MY_FREE(cur_path);
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
                   double *curmin, int *combi_min_idx, bool *skip_single,
                   int *active, int active_len, pthread_mutex_t *mutex) {
        
    
    
    double *a = MY_MALLOC(ncol * sizeof(double));
    CHECK_PTR(a, "a could not be allocated");
    double *b = MY_MALLOC(ncol * sizeof(double));
    CHECK_PTR(b, "b could not be allocated");
    double *c = MY_MALLOC(ncol * sizeof(double));
    CHECK_PTR(c, "c could not be allocated");
    double *d = MY_MALLOC(ncol * sizeof(double));
    CHECK_PTR(d, "d could not be allocated");
    double *gamma = MY_MALLOC(ncol * sizeof(double));
    CHECK_PTR(gamma, "gamma could not be allocated");

    double minval = DBL_MAX;
    int min_idx = -1;

    for (int j = 0; j < ncol; j++) {
        if(skip_single[j])
            continue;

        tot_node_explored++;

        a[j] = 0;
        b[j] = 0;
        c[j] = 0;
        d[j] = 0;
        for (int i = 0; i < nrow; i++) {
            if(cur_exprs[i + j*nrow]) {
                a[j] += w1p[i];
                b[j] += w1m[i];
                c[j] += w2p[i];
                d[j] += w2m[i];
            }
        }
                    
        gamma[j] = (num_const + a[j] - b[j]) / (denom_const + c[j] - d[j]);
        
        if(gamma[j] > EPSILON && gamma[j] < minval && !is_idx_in_set(cur_path, tot_depth, j, active, active_len)) {
            minval = gamma[j];
            min_idx = j;
        }
    }
//    printf("#%.0f | Single skipped: %ld\nskip_single[17375] = %d\n", tot_explored, single_skipped, skip_single[17375]);
    
//    int min_idx = get_pos_min_idx(gamma, ncol, cur_path, tot_depth, active, active_len);
    
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
        
        
        char *combi_exprs = MY_MALLOC(ncol*nrow);
        CHECK_PTR(combi_exprs, "combi_exprs could not be allocated");
        char *exprs_mask = MY_MALLOC(ncol*nrow);
        if(exprs_mask == NULL) {
            printf("exprs_mask could not be allocated: %d * %d = %d\ncurdepth: %d\n", ncol, nrow, ncol*nrow, cur_depth);
            exit(-1);
        }
        
        clock_t begin=clock();
        
        char null_block [nrow];
        bzero(null_block, nrow);
        
        char *x = MY_MALLOC(nrow);
        int *x_subset = MY_MALLOC(nrow * sizeof(int));
        int *ordered = MY_MALLOC(nrow * sizeof(int));
        int *ranking = MY_MALLOC(nrow * sizeof(int));

        
        for(int j = 0; j < ncol; j++) {
            if(skip_single[j]) {
                easy_pruned++;
                goto skip;
            }

#ifdef USE_HEURISTIC_1
            if((denom_const - d[j] >= 0
                && (((num_const - b[j]) / (denom_const + c[j]) >= *curmin) ||
                    ((num_const + a[j]) / (denom_const - d[j]) <= EPSILON)))
               || (denom_const + c[j] <= 0
                   && (((num_const + a[j]) / (denom_const - d[j]) >= *curmin) ||
                       ((num_const - b[j]) / (denom_const + c[j]) <= EPSILON)))) {
                       
                easy_pruned++;
                goto skip;
            }
#endif
            
            for(int i = 0; i < cur_depth; i++)
                if(cur_path[i] == j)
                    goto skip;

#ifdef USE_HEURISTIC_2

            if(denom_const - d[j] > 0 || denom_const + c[j] < 0) {
                
                bool pos_denom = (denom_const - d[j] > 0);
                double num, denom;
                
                int tot_subset = 0;
                                
                for (int i = 0; i < nrow; i++) {
                    
                    if(cur_exprs[i + j*nrow] > 0) {
                        
                        if(pos_denom) {
                            num = w1p[i] - w1m[i];
                            denom = w2p[i] - w2m[i];
                        }
                        else {
                            num = w1m[i] - w1p[i];
                            denom = w2m[i] - w2p[i];
                        }
                        
                        if((num < 0 && denom >= 0) || (num == 0 && denom > 0))
                            x[i] = 1;
                        else if((num > 0 && denom <= 0) || (num == 0 && denom < 0))
                            x[i] = 0;
                        else if(num == 0 && denom == 0)
                            x[i] = 0;
                        else {
                            x_subset[tot_subset++] = i;
                            
                            if(num < 0 && denom < 0) // both strictly negative -> equiv. to both positives, taking !x afterward
                                x[i] = -2;
                            else // both strictly positive
                                x[i] = -1;
                        }
                    }
                    else
                        x[i] = 0;
                }
                
                double optim_min = MAXFLOAT;
                
                if(tot_subset > 0)  {
                    
                    pthread_mutex_lock(mutex);
                    for (int i = 0; i < tot_subset; i++) {
                        int sub_i = x_subset[i];
                        if(pos_denom)
                            ordered[i] = tot_subset-1-i;
                        else
                            ordered[i] = i;
                        weight_ratios[i] = (w1m[sub_i] - w1p[sub_i])/(w2m[sub_i] - w2p[sub_i]);
                    }
                    qsort(ordered, tot_subset, sizeof(int), order_of_double_cmp);
                    pthread_mutex_unlock(mutex);

                    for (int i = 0; i < tot_subset; i++)
                        ranking[ordered[i]] = i;
                    
                    
                    
                    
                    double optim_min_n_found, optim_min_d_found;
//                    int save_k = 0;
                    
                    for (int i = -1; i <= tot_subset; i++) {
                        int sub_k = 0;
                        optim_min_n_found = num_const;
                        optim_min_d_found = denom_const;
                        
                        for (int k = 0; k < nrow; k++) {
                            if((x[k] > 0) || (x[k] == -1 && ranking[sub_k] <= i) || (x[k] == -2 && ranking[sub_k] > i)) {
                                optim_min_n_found += w1p[k] - w1m[k];
                                optim_min_d_found += w2p[k] - w2m[k];
                            }
                            
                            if(x[k] < 0)
                                sub_k++;
                        }
                        if(optim_min_n_found / optim_min_d_found < optim_min) {
                            optim_min = optim_min_n_found / optim_min_d_found;
//                            save_k = i;
                            if(optim_min < *curmin)
                                break;
                        }
                    }
                    
                }
                else {
                    double optim_min_n = num_const, optim_min_d = denom_const;
                    for (int i = 0; i < nrow; i++) {
                        if(x[i] > 0) {
                            optim_min_n += w1p[i] - w1m[i];
                            optim_min_d += w2p[i] - w2m[i];
                        }
                    }
                    optim_min = optim_min_n/optim_min_d;
                }

                
                if(optim_min >= *curmin) {
                    pruned++;
                    goto skip;
                }
                else
                    not_pruned++;
            }
#endif
            

            if(cur_depth > 0) { //doing breadth first for first level saves time
                for(int k = 0; k < ncol; k++)
                    memcpy(&exprs_mask[k*nrow], &cur_exprs[j*nrow], nrow);
                
                bool non_null = false;
                for (int i = 0; i < nrow*ncol / sizeof(uint64_t); i++) {
                    ((uint64_t *) combi_exprs)[i] = ((uint64_t *) single_exprs)[i] & ((uint64_t *) exprs_mask)[i];
                    non_null |= (((uint64_t *) combi_exprs)[i] != 0);
                }
                for (int i = nrow*ncol - (nrow*ncol % sizeof(uint64_t)); i < nrow*ncol; i++) {
                    combi_exprs[i] = single_exprs[i] & exprs_mask[i];
                    non_null |= (combi_exprs[i] != 0);
                }
                
                if (! non_null) {
                    easy_pruned++;
                    goto skip;
                }

                cur_path[cur_depth] = j;
                get_min_recur(cur_depth+1, tot_depth, cur_path, combi_exprs, single_exprs, nrow, ncol, w1p, w1m, w2p, w2m, num_const, denom_const, curmin, combi_min_idx, skip_single, active, active_len, mutex);
            }
            
            tot_non_leaf_explored++;

            continue;
            
            skip:
            if(cur_depth == 0)
                skip_single[j] = TRUE;
        }
        
        if (cur_depth == 0) {
            for (int j = 0; j < ncol; j++) {
                
                if(skip_single[j])
                    continue;
                
                for(int k = 0; k < ncol; k++)
                    memcpy(&exprs_mask[k*nrow], &cur_exprs[j*nrow], nrow);
                
                bool non_null = false;
                for (int i = 0; i < nrow*ncol / sizeof(uint64_t); i++) {
                    ((uint64_t *) combi_exprs)[i] = ((uint64_t *) single_exprs)[i] & ((uint64_t *) exprs_mask)[i];
                    non_null |= (((uint64_t *) combi_exprs)[i] != 0);
                }
                for (int i = nrow*ncol - (nrow*ncol % sizeof(uint64_t)); i < nrow*ncol; i++) {
                    combi_exprs[i] = single_exprs[i] & exprs_mask[i];
                    non_null |= (combi_exprs[i] != 0);
                }
                
                if (! non_null) {
                    easy_pruned++;
                    continue;
                }

                cur_path[cur_depth] = j;
                get_min_recur(cur_depth+1, tot_depth, cur_path, combi_exprs, single_exprs, nrow, ncol, w1p, w1m, w2p, w2m, num_const, denom_const, curmin, combi_min_idx, skip_single, active, active_len, mutex);

            }
        }
        
        
        cur_path[cur_depth] = -1;
        
        MY_FREE(ordered);
        MY_FREE(ranking);

        MY_FREE(x);
        MY_FREE(x_subset);

        MY_FREE(combi_exprs);
        MY_FREE(exprs_mask);
        
        clock_t end=clock();
        double diff = (end - begin)*1000/CLOCKS_PER_SEC;
        if(tot_non_leaf_explored > 0 && cur_depth < 1) {
            if(g_trace >= 2)
                printf("Explored %.0f non-terminal so far (total: %.2f s.).\n", tot_non_leaf_explored, diff/1000);
        }
    }
    
    MY_FREE(a);
    MY_FREE(b);
    MY_FREE(c);
    MY_FREE(d);
    MY_FREE(gamma);
    
    return;
}
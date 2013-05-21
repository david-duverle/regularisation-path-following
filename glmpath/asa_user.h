#include <limits.h>
#include <string.h>
//#include <stdlib.h>
#include <stdio.h>

#define INT int
#define INT_INF INT_MAX
#define INF DBL_MAX
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif
#ifndef NULL
#define NULL 0
#endif
/*============================================================================
   structure which is passed to the user's evaluation routines when
   either the objective function or gradient must be evaluated                */
typedef struct asa_objective_struct
{
    double      *x ; /* current iterate */
    double      *g ; /* user will store the gradient at x in g */
    INT          n ; /* problem dimension */
    int     *ifree ; /* if not NULL, contains indices of free variables */
    int      nfree ; /* if ifree not NULL, contains number of free variables */
    double *z;  /* user input */
    double *mz; /* user input */
} asa_objective ;

/*============================================================================
   user controlled parameters for gradient projection algorithm
               (default values in asa_default)                                */
typedef struct asa_parm_struct
{
    /* parameters values that the user may wish to modify */
/*----------------------------------------------------------------------------*/
    /* T => print final statistics
       F => no printout of statistics */
    int PrintFinal ;

    /* Level 0  = no printing, ... , Level 4 = maximum printing */
    int PrintLevel ;

    /* T => print parameters values
       F => do not display parmeter values */
    int PrintParms ;

    /* T => use approximate nonmonotone Armijo line search
       F => use ordinary nonmonotone Armijo line search, switch to
            approximate Armijo when |f_r-f| < AArmijoFac*|min (f_r, f_{max})| */
    int    AArmijo ;
    double AArmijoFac ;

    /* Stop Rules:
       T => ||proj_grad||_infty <= max(grad_tol,initial ||grad||_infty*StopFac)
       F => ||proj_grad||_infty <= grad_tol*(1 + |f_k|) */
    int    StopRule ;
    double StopFac ;

    /* T => estimated error in function value = eps*|min (f_r, f_{max}) |
       F => estimated error in function value = eps */
    int    PertRule ;
    double eps ;

    /* T => only use gradient projection algorithm
       F => let algorithm decide between grad_proj and cg_descent */
    int GradProjOnly;

    /* abort cbb when Armijo line search backtracks at least max_backstep
       times */
    int max_backsteps ;

    /* abort cbb after maxit_fac*n iterations in one pass through cbb */
    double maxit_fac ;

    /* abort cbb after totit_fac*n iterations in all passes through cbb */
    double totit_fac ;

    /* abort cbb iteration after maxfunc_fac*n function evaluations */
    double maxfunc_fac ;

    /* loj + pert_lo < xj < hij - pert_hi => xj free */
    double pert_lo ;
    double pert_hi ;

    /* search for non nan function value by shrinking search interval
       at most nshrink times */
    int nshrink ;

    /* factor by which interval shrinks when searching for non nan value */
    double nan_fac ;

/*============================================================================
       technical parameters which the user probably should not touch          */
    int                  L ; /* update fr if fmin was not improved after
                                L iterations */
    int                  m ; /* fmax = max (f_{k-i}, i = 0, 1, ...,
                                min (k, m-1) ) */
    int                  P ; /* update fr if P previous initial stepsize was
                                all accepted */
    int                 nm ; /* CBB cycle length */
    double           gamma ; /* criterion for reinitializing BB stepsize */
    double          gamma1 ; /* criterion for updating reference value fr */
    double          gamma2 ; /* criterion for updating reference value fr */
    double           delta ; /* Armijo line search parameter */
    double            lmin ; /* Lower bound for initial stepsize */
    double            lmax ; /* Upper bound for initial stepsize */
    double           parm1 ; /* used when attempting a quadratic
                                interpolation step */
    double           parm2 ; /* used when attempting a quadratic
                                interpolation step */
    double           parm3 ; /* criterion for reinitializing the BB stepsize */
    int              parm4 ; /* maximum previous BB steps used
                                when s^t y <= ZERO */
    double            tau1 ; /* if ginorm < tau1*pgnorm, continue gradient
                                projection steps  */
    double      tau1_decay ; /* decay factor for tau1 */
    double            tau2 ; /* ginorm < tau2*pgnorm implies subproblem
                                solved in cgdescent */
    double      tau2_decay ; /* decay factor for tau2 */
    double         pgdecay ; /* criterion for checking undecided index set */
    double    armijo_decay ; /* decay factor in Armijo line search */
    double         armijo0 ; /* criterion for quadratic interpolation in
                                cbb line search */
    double         armijo1 ; /* criterion for quadratic interpolation in
                                cbb line search */
} asa_parm ;

/*============================================================================
   user controlled parameters for the conjugate gradient algorithm
               (default values in asa_cg_default)                             */
typedef struct asacg_parm_struct
{
    /* parameters values that the user may wish to modify */
/*----------------------------------------------------------------------------*/
    /* Level 0  = no printing), ... , Level 4 = maximum printing */
    int PrintLevel ;

    /* T => print parameters values
       F => do not display parmeter values */
    int PrintParms ;

    /* T => use approximate Wolfe line search
       F => use ordinary Wolfe line search, switch to approximate Wolfe when
                |f_k+1-f_k| < AWolfeFac*C_k, C_k = average size of cost  */
    int    AWolfe ;
    double AWolfeFac ;

    /* factor in [0, 1] used to compute average cost magnitude C_k as follows:
       Q_k = 1 + (Qdecay)Q_k-1, Q_0 = 0,  C_k = C_k-1 + (|f_k| - C_k-1)/Q_k */
    double Qdecay ;

    /* T => estimated error in function value is eps*Ck,
       F => estimated error in function value is eps */
    int    PertRule ;
    double eps ;

    /* T => attempt quadratic interpolation in line search when
                |f_k+1 - f_k|/f_k <= QuadCutOff
       F => no quadratic interpolation step */
    int    QuadStep ;
    double QuadCutOff ;

    /* T => check that f_k+1 - f_k <= debugtol*C_k
       F => no checking of function values */
    int    debug ;
    double debugtol ;

    /* if step is nonzero, it is the initial step of the initial line search */
    double step ;

    /* abort cg after maxit_fac*n iterations in one pass */
    double maxit_fac ;

    /* abort cg after totit_fac*n iterations in all passes */
    double totit_fac ;

    /* maximum number of times the bracketing interval grows or shrinks
       in the line search is nexpand */
    int nexpand ;

   /* maximum number of secant iterations in line search is nsecant */
    int nsecant ;

    /* conjugate gradient method restarts after (n*restart_fac) iterations */
    double restart_fac ;

    /* stop when -alpha*dphi0 (estimated change in function value) <= feps*|f|*/
    double feps ;

    /* after encountering nan, growth factor when searching for
       a bracketing interval */
    double nan_rho ;

/*============================================================================
       technical parameters which the user probably should not touch          */
    double           delta ; /* Wolfe line search parameter */
    double           sigma ; /* Wolfe line search parameter */
    double           gamma ; /* decay factor for bracket interval width */
    double             rho ; /* growth factor when searching for initial
                                bracketing interval */
    double             eta ; /* lower bound for the conjugate gradient update
                                parameter beta_k is eta*||d||_2 */
    double            psi0 ; /* factor used in starting guess for iteration 1 */
    double            psi1 ; /* in performing a QuadStep, we evaluate the
                                function at psi1*previous step */
    double            psi2 ; /* when starting a new cg iteration, our initial
                                guess for the line search stepsize is
                                psi2*previous step */
} asacg_parm ;

typedef struct asa_stat_struct /* statistics returned to user */
{
    double               f ; /*function value at solution */
    double          pgnorm ; /* ||Proj (x_k - g_k) - x_k||_infty */
    INT            cbbiter ; /* total cbb iterations */
    INT            cbbfunc ; /* total cbb function evaluations */
    INT            cbbgrad ; /* total cbb gradient evaluations */
    INT             cgiter ; /* total cg iterations */
    INT             cgfunc ; /* total cg function evaluations */
    INT             cggrad ; /* total cg gradient evaluations */
} asa_stat ;

/* prototypes */

int asa_cg /*  return:
                      -2 (function value became nan in cg)
                      -1 (starting function value is nan in cg)
                       0 (convergence tolerance satisfied)
                       1 (change in func <= feps*|f| in cg)
                       2 (cg iterations in all passes or
                          in one pass exceeded their limit)
                       3 (slope always negative in line search in cg)
                       4 (number secant iterations exceed nsecant in cg)
                       5 (search direction not a descent direction in cg)
                       6 (line search fails in initial interval in cg)
                       7 (line search fails during bisection in cg)
                       8 (line search fails during interval update in cg)
                       9 (debugger is on and the function value increases in cg)
                      10 (out of memory)
                      11 (cbb iterations in all passes or
                          in one pass exceeded their limit)
                      12 (line search failed in cbb iteration)
                      13 (search direction in cbb is not descent direction)
                      14 (function value became nan in cbb) */
(
    double            *x, /* input: starting guess, output: the solution */
    double           *lo, /* lower bounds */
    double           *hi, /* upper bounds */
    INT                n, /* problem dimension */
    asa_stat       *Stat, /* structure with statistics (can be NULL) */
    asacg_parm    *CParm, /* user parameters, NULL = use default parameters */
    asa_parm      *AParm, /* user parameters, NULL = use default parameters */
    double      grad_tol, /* |Proj (x_k - g_k) - x_k|_inf <= grad_tol */
    double   (*value) (asa_objective *), /* evaluate objective function */
    void      (*grad) (asa_objective *), /* evaluate objective gradient */
    double (*valgrad) (asa_objective *), /* function and gradient
                                            NULL = use value & grad routines */
    double        *Work,  /* either work array of size 7n + memory (m) or NULL */
    double *z,
    double *mz
) ;

void asa_default /* set default parameter values for asa */
(
    asa_parm   *Parm
) ;

void asa_cg_default /* set parameter values for cg_descent */
(
    asacg_parm   *Parm
) ;

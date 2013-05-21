#include <R.h>
#include <math.h>
#include "asa_user.h"

double glmpath_value
(
    asa_objective *asa
) ;

void glmpath_grad
(
    asa_objective *asa
) ;

double glmpath_valgrad
(
    asa_objective *asa
) ;

void solve_glmpath (int *n, double *x,
                    double *lo, double *hi,
                    int *solved,
                    double *z, double *mz)
{
  asacg_parm cgParm ;
  asa_parm asaParm ;

  asa_cg_default (&cgParm) ;
  asa_default (&asaParm) ;
  cgParm.PrintLevel = 0 ;
  asaParm.PrintLevel = 0 ;

  *solved = asa_cg (x, lo, hi, *n, NULL, &cgParm, &asaParm, 1.e-8,
                    glmpath_value, glmpath_grad, glmpath_valgrad,
                    NULL, z, mz) ;
  return ;
}

double glmpath_value
(
    asa_objective *asa
)
{
  // z: data_x(m*p), data_y(m), w(m), offset(m), penalty(p)
  // mz: m(nobs), distr, lam1, lam2
  int n;
  double f, *x, *z, *mz;
  int p, m, distr, iy, iw, ioffset, ipenalty, i, j;
  double *b, *y, *w, *eta, *mu;
  double lamone, lamtwo, loglik, normone, normtwo;

  n = asa->n;
  x = asa->x;
  z = asa->z;
  mz = asa->mz;
  f = 0.0;

  p = n / 2;
  m = mz[0];
  distr = mz[1];
  lamone = mz[2];
  lamtwo = mz[3];

  b  = (double *) malloc (p * sizeof (double));
  y  = (double *) malloc (m * sizeof (double));
  w  = (double *) malloc (m * sizeof (double));
  eta  = (double *) malloc (m * sizeof (double));
  mu  = (double *) malloc (m * sizeof (double));

  // b = x(+) - x(-); b[j] = x[j] - x[p + j]
  for (j = 0; j < p; ++j) {
    b[j] = x[j] - x[p + j];
  }

  iy = p * m;
  iw = iy + m;
  ioffset = iw + m;
  ipenalty = ioffset + m;

  for (i = 0; i < m; ++i) {
    y[i] = z[iy + i];
    w[i] = z[iw + i];
    eta[i] = z[ioffset + i];
    for (j = 0; j < p; ++j) {
      eta[i] += (b[j] * z[j * m + i]);
    }
  }

  loglik = 0.0;
  for (i = 0; i < m; ++i) {
    if (distr == 0) {
      // normal distribution
      mu[i] = eta[i];
      loglik -= (w[i] * pow((y[i] - mu[i]), 2) / 2.0);
    } else if (distr == 1) {
      // binomial distribution
      mu[i] = 1.0 / (1.0 + exp(-1.0 * eta[i]));
      loglik += (w[i] * (y[i] * eta[i] -
                         log(1.0 + exp(eta[i]))));
    } else if (distr == 2) {
      // poisson distribution
      mu[i] = exp(eta[i]);
      loglik += (w[i] * (y[i] * eta[i] - mu[i]));
    }
  }

  normone = 0.0;
  normtwo = 0.0;
  for (j = 0; j < p; ++j) {
    if (z[ipenalty + j] == 1.0) {
      normone += fabs(b[j]);
    }
    if (j > 0) {
      normtwo += pow(b[j], 2);
    }
  }

  f = (-1.0 * loglik +
       lamone * normone + lamtwo * normtwo / 2.0);

  free(b);
  free(y);
  free(w);
  free(eta);
  free(mu);
  return (f);
}

void glmpath_grad
(
    asa_objective *asa
)
{
  int n;
  double *g, *x, *z, *mz;
  int p, m, distr, iy, iw, ioffset, ipenalty, i, j;
  double *b, *y, *w, *eta, *mu, *resid;
  double lamone, lamtwo;

  n = asa->n;
  x = asa->x;
  g = asa->g;
  z = asa->z;
  mz = asa->mz;

  p = n / 2;
  m = mz[0];
  distr = mz[1];
  lamone = mz[2];
  lamtwo = mz[3];

  b  = (double *) malloc (p * sizeof (double));
  y  = (double *) malloc (m * sizeof (double));
  w  = (double *) malloc (m * sizeof (double));
  eta  = (double *) malloc (m * sizeof (double));
  mu  = (double *) malloc (m * sizeof (double));
  resid = (double *) malloc (m * sizeof(double));

  // b = x(+) - x(-); b[j] = x[j] - x[p + j]
  for (j = 0; j < p; ++j) {
    b[j] = x[j] - x[p + j];
  }

  iy = p * m;
  iw = iy + m;
  ioffset = iw + m;
  ipenalty = ioffset + m;

  for (i = 0; i < m; ++i) {
    y[i] = z[iy + i];
    w[i] = z[iw + i];
    eta[i] = z[ioffset + i];
    for (j = 0; j < p; ++j) {
      eta[i] += (b[j] * (z[j * m + i]));
    }
  }

  for (i = 0; i < m; ++i) {
    if (distr == 0) {
      // normal distribution
      mu[i] = eta[i];
    } else if (distr == 1) {
      // binomial distribution
      mu[i] = 1.0 / (1.0 + exp(-1.0 * eta[i]));
    } else if (distr == 2) {
      // poisson distribution
      mu[i] = exp(eta[i]);
    }
    resid[i] = w[i] * (y[i] - mu[i]);
  }

  for (j = 0; j < p; ++j) {
    g[j] = 0.0;
    for (i = 0; i < m; ++i) {
      g[j] -= (z[j * m + i] * resid[i]);
    }
    g[p + j] = -1.0 * g[j];
    if (z[ipenalty + j] == 1.0) {
      g[j] += lamone;
      g[p + j] += lamone;
    }
    if (j > 0) {
      g[j] += (lamtwo * b[j]);
      g[p + j] -= (lamtwo * b[j]);
    }
  }

  free(b);
  free(y);
  free(w);
  free(eta);
  free(mu);
  free(resid);
  return ;
}

double glmpath_valgrad
(
    asa_objective *asa
)
{
  // z: data_x(m*p), data_y(m), w(m), offset(m), penalty(p)
  // mz: m(nobs), distr, lam1, lam2
  int n;
  double f, *g, *x, *z, *mz;
  int p, m, distr, iy, iw, ioffset, ipenalty, i, j;
  double *b, *y, *w, *eta, *mu, *resid;
  double lamone, lamtwo, loglik, normone, normtwo;

  n = asa->n;
  x = asa->x;
  g = asa->g;
  z = asa->z;
  mz = asa->mz;
  f = 0.0;

  p = n / 2;
  m = mz[0];
  distr = mz[1];
  lamone = mz[2];
  lamtwo = mz[3];

  b  = (double *) malloc (p * sizeof (double));
  y  = (double *) malloc (m * sizeof (double));
  w  = (double *) malloc (m * sizeof (double));
  eta  = (double *) malloc (m * sizeof (double));
  mu  = (double *) malloc (m * sizeof (double));
  resid = (double *) malloc (m * sizeof(double));

  // b = x(+) - x(-); b[j] = x[j] - x[p + j]
  for (j = 0; j < p; ++j) {
    b[j] = x[j] - x[p + j];
  }

  iy = p * m;
  iw = iy + m;
  ioffset = iw + m;
  ipenalty = ioffset + m;

  for (i = 0; i < m; ++i) {
    y[i] = z[iy + i];
    w[i] = z[iw + i];
    eta[i] = z[ioffset + i];
    for (j = 0; j < p; ++j) {
      eta[i] += (b[j] * z[j * m + i]);
    }
  }

  loglik = 0.0;
  for (i = 0; i < m; ++i) {
    if (distr == 0) {
      // normal distribution
      mu[i] = eta[i];
      loglik -= (w[i] * pow((y[i] - mu[i]), 2) / 2.0);
    } else if (distr == 1) {
      // binomial distribution
      mu[i] = 1.0 / (1.0 + exp(-1.0 * eta[i]));
      loglik += (w[i] * (y[i] * eta[i] -
                         log(1.0 + exp(eta[i]))));
    } else if (distr == 2) {
      // poisson distribution
      mu[i] = exp(eta[i]);
      loglik += (w[i] * (y[i] * eta[i] - mu[i]));
    }
    resid[i] = w[i] * (y[i] - mu[i]);
  }

  for (j = 0; j < p; ++j) {
    g[j] = 0.0;
    for (i = 0; i < m; ++i) {
      g[j] -= (z[j * m + i] * resid[i]);
    }
    g[p + j] = -1.0 * g[j];
    if (z[ipenalty + j] == 1.0) {
      g[j] += lamone;
      g[p + j] += lamone;
    }
    if (j > 0) {
      g[j] += (lamtwo * b[j]);
      g[p + j] -= (lamtwo * b[j]);
    }
  }

  normone = 0.0;
  normtwo = 0.0;
  for (j = 0; j < p; ++j) {
    if (z[ipenalty + j] == 1.0) {
      normone += fabs(b[j]);
    }
    if (j > 0) {
      normtwo += pow(b[j], 2);
    }
  }

  f = (-1.0 * loglik +
       lamone * normone + lamtwo * normtwo / 2.0);

  free(b);
  free(y);
  free(w);
  free(eta);
  free(mu);
  free(resid);
  return (f);
}

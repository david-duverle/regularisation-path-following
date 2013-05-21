#include <R.h>
#include <math.h>
#include "asa_user.h"

double coxpath_value
(
    asa_objective *asa
) ;

void coxpath_grad
(
    asa_objective *asa
) ;

double coxpath_valgrad
(
    asa_objective *asa
) ;

void solve_coxpath (int *n, double *x,
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
                    coxpath_value, coxpath_grad, coxpath_valgrad,
                    NULL, z, mz) ;
  return ;
}

double coxpath_value
(
    asa_objective *asa
)
{
  // z: data_x(m*p), d(m), rept(m), penalty(p), <eta>(m), <wsum>(m)
  // mz: m(nobs), method, lam1, lam2, <lp>
  int n;
  double f, *x, *z, *mz;
  int p, m, method;
  int id, irept, ipenalty, ieta, iwsum, ilp, i, j, k0, k, l;
  int *d, *rept;
  double *b, *eta, *mu;
  double lamone, lamtwo, normone, normtwo, wsum, dwsum;

  n = asa->n;
  x = asa->x;
  z = asa->z;
  mz = asa->mz;
  f = 0.0;

  p = n / 2;
  m = mz[0];
  method = mz[1];
  lamone = mz[2];
  lamtwo = mz[3];

  b  = (double *) malloc (p * sizeof (double));
  d  = (int *) malloc (m * sizeof (int));
  rept  = (int *) malloc (m * sizeof (int));
  eta  = (double *) malloc (m * sizeof (double));
  mu  = (double *) malloc (m * sizeof (double));

  // b = x(+) - x(-); b[j] = x[j] - x[p + j]
  for (j = 0; j < p; ++j) {
    b[j] = x[j] - x[p + j];
  }

  id = p * m;
  irept = id + m;
  ipenalty = irept + m;
  ieta = ipenalty + p;
  iwsum = ieta + m;
  ilp = 4;

  for (i = 0; i < m; ++i) {
    d[i] = z[id + i];
    rept[i] = z[irept + i];
    mu[i] = 0.0;
    for (j = 0; j < p; ++j) {
      mu[i] += (b[j] * z[j * m + i]);
    }
    eta[i] = exp(mu[i]);
    z[ieta + i] = eta[i];
  }

  k = 0;
  for (i = 0; i < m; ++i) {
    if (d[i] == 1) {
      if (method == 1) {
        if (rept[i] != 0 && k == 0) {
          wsum = 0.0;
          k = rept[i] - 1;
          for (l = 0; l <= (i + k); ++l) {
            wsum += eta[l];
          }
        } else if (k > 0) {
          k -= 1;
        }
      } else if (method == 2) {
        if (rept[i] != 0 && k == 0) {
          k0 = rept[i];
          k = rept[i] - 1;
          wsum = 0.0;
          dwsum = 0.0;
          for (l = 0; l <= (i + k); ++l) {
            wsum += eta[l];
            if (k > 0 && l >= i) {
              dwsum += eta[l];
            }
          }
        } else if (k > 0) {
          wsum -= (dwsum / k0);
          k -= 1;
        }
      }
      z[iwsum + i] = wsum;
      f += (-1.0 * mu[i] + log(wsum));
    }
  }
  mz[ilp] = -1.0 * f;

  normone = 0.0;
  normtwo = 0.0;
  for (j = 0; j < p; ++j) {
    if (z[ipenalty + j] == 1.0) {
      normone += fabs(b[j]);
    }
    normtwo += pow(b[j], 2);
  }

  f += (lamone * normone + lamtwo * normtwo / 2.0);

  free(b);
  free(d);
  free(rept);
  free(eta);
  free(mu);
  return (f);
}

void coxpath_grad
(
    asa_objective *asa
)
{
  // z: data_x(m*p), d(m), rept(m), penalty(p), <eta>(m), <wsum>(m)
  // mz: m(nobs), method, lam1, lam2, <lp>
  int n;
  double *g, *x, *z, *mz;
  int p, m, method;
  int iy, id, irept, ipenalty, i, j, k0, k, l;
  int *d, *rept;
  double *b, *eta, *mu, *wx, *dwx;
  double lamone, lamtwo, wsum, dwsum;

  n = asa->n;
  x = asa->x;
  g = asa->g;
  z = asa->z;
  mz = asa->mz;

  p = n / 2;
  m = mz[0];
  method = mz[1];
  lamone = mz[2];
  lamtwo = mz[3];

  b  = (double *) malloc (p * sizeof (double));
  d  = (int *) malloc (m * sizeof (int));
  rept  = (int *) malloc (m * sizeof (int));
  eta  = (double *) malloc (m * sizeof (double));
  mu  = (double *) malloc (m * sizeof (double));
  wx  = (double *) malloc (p * sizeof (double));
  dwx  = (double *) malloc (p * sizeof (double));

  // b = x(+) - x(-); b[j] = x[j] - x[p + j]
  for (j = 0; j < p; ++j) {
    b[j] = x[j] - x[p + j];
  }

  for (j = 0; j < n; ++j) {
    g[j] = 0.0;
  }

  id = p * m;
  irept = id + m;
  ipenalty = irept + m;

  for (i = 0; i < m; ++i) {
    d[i] = z[id + i];
    rept[i] = z[irept + i];
    mu[i] = 0.0;
    for (j = 0; j < p; ++j) {
      mu[i] += (b[j] * z[j * m + i]);
    }
    eta[i] = exp(mu[i]);
  }

  k = 0;
  for (i = 0; i < m; ++i) {
    if (d[i] == 1) {
      if (method == 1) {
        if (rept[i] != 0 && k == 0) {
          for (j = 0; j < p; ++j) {
            wx[j] = 0.0;
          }
          wsum = 0.0;
          k = rept[i] - 1;
          for (l = 0; l <= (i + k); ++l) {
            for (j = 0; j < p; ++j) {
              wx[j] += (z[j * m + l] * eta[l]);
            }
            wsum += eta[l];
          }
          for (j = 0; j < p; ++j) {
            wx[j] /= wsum;
          }
        } else if (k > 0) {
          k -= 1;
        }
      } else if (method == 2) {
        if (rept[i] != 0 && k == 0) {
          for (j = 0; j < p; ++j) {
            wx[j] = 0.0;
            dwx[j] = 0.0;
          }
          wsum = 0.0;
          dwsum = 0.0;
          k0 = rept[i];
          k = rept[i] - 1;
          for (l = 0; l <= (i + k); ++l) {
            for (j = 0; j < p; ++j) {
              wx[j] += (z[j * m + l] * eta[l]);
              if (k > 0 && l >= i) {
                dwx[j] += (z[j * m + l] * eta[l]);
              }
            }
            wsum += eta[l];
            if (k > 0 && l >= i) {
              dwsum += eta[l];
            }
          }
          for (j = 0; j < p; ++j) {
            wx[j] /= wsum;
          }
        } else if (k > 0) {
          for (j = 0; j < p; ++j) {
            wx[j] *= wsum;
            wx[j] -= (dwx[j] / k0);
          }
          wsum -= (dwsum / k0);
          for (j = 0; j < p; ++j) {
            wx[j] /= wsum;
          }
          k -= 1;
        }
      }
      for (j = 0; j < p; j++) {
        g[j] += (-1.0 * z[j * m + i] + wx[j]);
      }
    }
  }

  for (j = 0; j < p; ++j) {
    g[p + j] = -1.0 * g[j];
    if (z[ipenalty + j] == 1.0) {
      g[j] += lamone;
      g[p + j] += lamone;
    }
    g[j] += (lamtwo * b[j]);
    g[p + j] -= (lamtwo * b[j]);
  }

  free(b);
  free(d);
  free(rept);
  free(eta);
  free(mu);
  free(wx);
  free(dwx);
  return ;
}

double coxpath_valgrad
(
    asa_objective *asa
)
{
  // z: data_x(m*p), d(m), rept(m), penalty(p), <eta>(m), <wsum>(m)
  // mz: m(nobs), method, lam1, lam2, <lp>
  int n;
  double f, *g, *x, *z, *mz;
  int p, m, method;
  int iy, id, irept, ipenalty, ieta, iwsum, ilp, i, j, k0, k, l;
  int *d, *rept;
  double *b, *eta, *mu, *wx, *dwx;
  double lamone, lamtwo, normone, normtwo, wsum, dwsum;

  n = asa->n;
  x = asa->x;
  g = asa->g;
  z = asa->z;
  mz = asa->mz;
  f = 0.0;

  p = n / 2;
  m = mz[0];
  method = mz[1];
  lamone = mz[2];
  lamtwo = mz[3];

  b  = (double *) malloc (p * sizeof (double));
  d  = (int *) malloc (m * sizeof (int));
  rept  = (int *) malloc (m * sizeof (int));
  eta  = (double *) malloc (m * sizeof (double));
  mu  = (double *) malloc (m * sizeof (double));
  wx  = (double *) malloc (p * sizeof (double));
  dwx  = (double *) malloc (p * sizeof (double));

  // b = x(+) - x(-); b[j] = x[j] - x[p + j]
  for (j = 0; j < p; ++j) {
    b[j] = x[j] - x[p + j];
  }

  for (j = 0; j < n; ++j) {
    g[j] = 0.0;
  }

  id = p * m;
  irept = id + m;
  ipenalty = irept + m;
  ieta = ipenalty + p;
  iwsum = ieta + m;
  ilp = 4;

  for (i = 0; i < m; ++i) {
    d[i] = z[id + i];
    rept[i] = z[irept + i];
    mu[i] = 0.0;
    for (j = 0; j < p; ++j) {
      mu[i] += (b[j] * z[j * m + i]);
    }
    eta[i] = exp(mu[i]);
    z[ieta + i] = eta[i];
  }

  k = 0;
  for (i = 0; i < m; ++i) {
    if (d[i] == 1) {
      if (method == 1) {
        if (rept[i] != 0 && k == 0) {
          for (j = 0; j < p; ++j) {
            wx[j] = 0.0;
          }
          wsum = 0.0;
          k = rept[i] - 1;
          for (l = 0; l <= (i + k); ++l) {
            for (j = 0; j < p; ++j) {
              wx[j] += (z[j * m + l] * eta[l]);
            }
            wsum += eta[l];
          }
          for (j = 0; j < p; ++j) {
            wx[j] /= wsum;
          }
        } else if (k > 0) {
          k -= 1;
        }
      } else if (method == 2) {
        if (rept[i] != 0 && k == 0) {
          for (j = 0; j < p; ++j) {
            wx[j] = 0.0;
            dwx[j] = 0.0;
          }
          wsum = 0.0;
          dwsum = 0.0;
          k0 = rept[i];
          k = rept[i] - 1;
          for (l = 0; l <= (i + k); ++l) {
            for (j = 0; j < p; ++j) {
              wx[j] += (z[j * m + l] * eta[l]);
              if (k > 0 && l >= i) {
                dwx[j] += (z[j * m + l] * eta[l]);
              }
            }
            wsum += eta[l];
            if (k > 0 && l >= i) {
              dwsum += eta[l];
            }
          }
          for (j = 0; j < p; ++j) {
            wx[j] /= wsum;
          }
        } else if (k > 0) {
          for (j = 0; j < p; ++j) {
            wx[j] *= wsum;
            wx[j] -= (dwx[j] / k0);
          }
          wsum -= (dwsum / k0);
          for (j = 0; j < p; ++j) {
            wx[j] /= wsum;
          }
          k -= 1;
        }
      }
      z[iwsum + i] = wsum;
      for (j = 0; j < p; j++) {
        g[j] += (-1.0 * z[j * m + i] + wx[j]);
      }
      f += (-1.0 * mu[i] + log(wsum));
    }
  }
  mz[ilp] = -1.0 * f;

  for (j = 0; j < p; ++j) {
    g[p + j] = -1.0 * g[j];
    if (z[ipenalty + j] == 1.0) {
      g[j] += lamone;
      g[p + j] += lamone;
    }
    g[j] += (lamtwo * b[j]);
    g[p + j] -= (lamtwo * b[j]);
  }

  normone = 0.0;
  normtwo = 0.0;
  for (j = 0; j < p; ++j) {
    if (z[ipenalty + j] == 1.0) {
      normone += fabs(b[j]);
    }
    normtwo += pow(b[j], 2);
  }

  f += (lamone * normone + lamtwo * normtwo / 2.0);

  free(b);
  free(d);
  free(rept);
  free(eta);
  free(mu);
  free(wx);
  free(dwx);
  return (f);
}

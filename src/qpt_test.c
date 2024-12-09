#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

void _Gram(double *x, int *param, double *xbyx, double **xTx)
{
    // input: 
    // x in R^{n*p}, saved by row.
    // param = c(n, p).
    // xbyx is the upper triangle of x^T x without diagonal, saved by row.
    // xTx is the matrix x^T x, 
    // only saves the upper triangle without diagonal, saved by row.
    unsigned int n, p, i, j, k, c;
    n = param[0];
    p = param[1];
    c = 1;
    double tmp, *xi, *xj;

    for (i = 0; i < n-1; ++i) {
        xi = x + i*p; // point to i-th sample
        for (j = i+1; j < n; ++j) {
            xj = x + j*p; // point to j-th sample
            for (tmp = 0.0, k = 0; k < p; ++k) 
                tmp += xi[k]*xj[k];
            xbyx[c] = tmp;
            c++;
        }
    }

    c = 0;
    for (i = 0; i < n-1; ++i) {
        xTx[i] = xbyx + c;
        c += n-2-i;
    }
}

double _qpt_test(double **xTx, double *Ytau, int n)
{   
    unsigned int i, j, k, ell;
    double Test = 0.0;
    double Test0 = 0.0;
    double tmp = (n-2.0)*(n-3.0);
    for (i = 0; i < n-3; ++i) {
        for (j = i+1; j < n-2; ++j) {
            Test0 = 0.0;
            for (k = j+1; k < n-1; ++k) {
                for (ell = k+1; ell < n; ++ell) {
                    Test0 += (Ytau[i] - Ytau[j]) * (Ytau[k] - Ytau[ell]) * (xTx[i][k] - xTx[j][k] - xTx[i][ell] + xTx[j][ell]);
                    Test0 += (Ytau[i] - Ytau[k]) * (Ytau[j] - Ytau[ell]) * (xTx[i][j] - xTx[j][k] - xTx[i][ell] + xTx[k][ell]);
                    Test0 += (Ytau[i] - Ytau[ell]) * (Ytau[j] - Ytau[k]) * (xTx[i][j] - xTx[j][ell] - xTx[i][k] + xTx[k][ell]);
                }
            }
            Test += Test0 * 8.0 / tmp;
        }
    }
    Test /= n*(n-1.0);

    return Test;
}


double _tr_sigma2_hat(double **xTx, int n)
{
    int i, j, k, l;
    double trsigma2, tmp, y1n, y2n, y3n;
    y1n = y2n = y3n = 0.0;

    for (i = 0; i < n-1; ++i) {
        for (j = i+1; j < n; ++j) {
            y1n += xTx[i][j] * xTx[i][j];
            for (k = j+1; k < n; ++k) {
                tmp = (xTx[i][j] + xTx[i][k])*xTx[j][k] + xTx[i][j]*xTx[i][k];
                y2n += tmp/(n-2);
                for (l = k+1; l < n; ++l) {
                    tmp = xTx[i][j]*xTx[k][l] + xTx[i][k]*xTx[j][l] + xTx[i][l]*xTx[j][k];
                    y3n += tmp*4/(n-2)/(n-3);
                }
            }
        }
    }

    trsigma2 = (y1n - y2n*2 + y3n)/(n*(n-1)/2);
    return trsigma2;
}


SEXP _QPT_Test(SEXP _X, SEXP _YTAU, SEXP _Param, SEXP _TAU)
{
    int n;
    n = INTEGER(_Param)[0];
    double *ytau;
    ytau = REAL(_YTAU);
    double tau;
    tau = REAL(_TAU)[0];

    double *xbyx;
    double **xTx;
    xbyx     =  (double*)malloc(sizeof(double) * (n*n-n+2)/2);
    xTx      = (double**)malloc(sizeof(double*)*(n-1));

    SEXP _Test, _Tr_sigma2;
    PROTECT(_Test = allocVector(REALSXP, 1));
    PROTECT(_Tr_sigma2 = allocVector(REALSXP, 1));

    double test, tr_sg2, test_sd;

    _Gram(REAL(_X), INTEGER(_Param), xbyx, xTx);
    test = _qpt_test(xTx, ytau, n);
    tr_sg2 = _tr_sigma2_hat(xTx, n);
    test_sd = 4.0*tau*(1.0-tau)*sqrt(2.0*tr_sg2) / n;
    test /= test_sd;

    REAL(_Test)[0] = test;
    REAL(_Tr_sigma2)[0] = tr_sg2;

    SEXP Result, R_names;
    char *names[2] = {"test", "tr_sigma"};
    PROTECT(Result    = allocVector(VECSXP,  2));
    PROTECT(R_names   = allocVector(STRSXP,  2));
    SET_STRING_ELT(R_names, 0,  mkChar(names[0]));
    SET_STRING_ELT(R_names, 1,  mkChar(names[1]));
    SET_VECTOR_ELT(Result, 0, _Test);
    SET_VECTOR_ELT(Result, 1, _Tr_sigma2);
    setAttrib(Result, R_NamesSymbol, R_names); 

    free(xbyx);
    free(xTx);
    UNPROTECT(4);
    return Result;
}



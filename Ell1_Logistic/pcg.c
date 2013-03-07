/** \file   pcg.c
 *  \brief  Source file for a preconditioned conjugate gradient method
 *          implementation.
 *
 *  This PCG implementation is based on the algorithm in the book
 *  Numerical Optimization by Nocedal and Write.
 */

#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef Rpackage
    #include "blas.h"
#else
    #include <R_ext/BLAS.h>
#endif


#include "def.h"
#include "pcg.h"

/** \brief Initialize pcg_status_t structure.
 *
 *  @param  pcgstat pointer to a pcg_status_t structure.
 */
void init_pcg_status(pcg_status_t *pcgstat)
{
    pcgstat->flag   = PCG_INIT;
    pcgstat->iter   = 0;
    pcgstat->relres = 0.0;
}

/** \brief Solves Ax = b via PCG.
 *
 *  PCG implements a preconditioned conjugate gradient method,
 *  which solves a positive definite linear system Ax = b.
 *  Here, A is a (m x n) matrix, b is an m-vector, and 
 *  variable x is an n-vector.
 *
 *  @param  x       double array of size n;
 *                  This parameter is used for both taking an initial value
 *                  and returning the result.
 *  @param  pcgstat pointer to a structure of PCG status.
 *  @param  afun    pointer to a function that computes A*x.
 *  @param  adata   pointer to a data structure for afun;
 *  @param  mfun    pointer to a function that computes M^{-1}r.
 *  @param  mdata   pointer to a data structure for mfun;
 *  @param  b       pointer to a n-vector
 *  @param  tol     integer; terminates algorithm when ||r|| < ||b||*tol.
 *  @param  maxiter integer; terminates algorithm when iteration 
 *                  exceeds maxiter.
 *  @param  n       integer; dimension of the problem.
 *
 *  Call back functions:
 *  - afun(double *out, const double *in, const void *adata):
 *  -- afun is a function which computes out = A*in.
 *  - mfun(double *out, const double *in, const void *mdata)
 *  -- mfun is a function which computes out = M^(-1)*in.
 */
int pcg(double *x, pcg_status_t *pcgstat,
        void (*afun)(double*, const double*, const void*), const void* adata,
        void (*mfun)(double*, const double*, const void*), const void* mdata,
        double *b, double tol, double maxiter, int n)
{
    int    k;       /* pcg iteration                */
    double alpha;   /* x_{k+1} = x_k + alpha_k*p_k  */
    double beta;    /* p_k = -r_k * beta_k*p_{k-1}  */
    double bnorm;   /* ||b||                        */
    double rnorm;   /* ||r||                        */

    double *r;      /* Residual                     */
    double *y;      /* Temporary for Ap,M^{-1}r     */
    double *p;      /* Conjugate direction          */
    double ry;      /* r^T*y                        */
    double ry2;

    int    mink;
    int    nbytes;  /* number of bytes of a vector  */
    nbytes  = n*sizeof(double);

    r = malloc(nbytes);
    y = malloc(nbytes);
    p = malloc(nbytes);

    /* r := Ax - b */
    afun(r, x, adata);
    F77_CALL(daxpy)(&n, &dminusone, b, &ione, r, &ione);

    /* y := M^{-1}r */
    mfun(y, r, mdata);

    /* p := -y */
    memcpy(p, y, nbytes);
    F77_CALL(dscal)(&n, &dminusone, p, &ione);

    ry    = F77_CALL(ddot)(&n, r, &ione, y, &ione);
    rnorm = F77_CALL(dnrm2)(&n, r, &ione);
    bnorm = F77_CALL(dnrm2)(&n, b, &ione);

    mink = min(n/2,10);
    
    for (k = 0; k < mink || (k<maxiter && (rnorm=F77_CALL(dnrm2)(&n, r, &ione))>bnorm*tol) ; k++)
    {
        afun(y, p, adata);
        alpha = ry / F77_CALL(ddot)(&n, p, &ione, y, &ione); /* alpha := (r^Ty)/(p^TAp)  */

        F77_CALL(daxpy)(&n, &alpha, p, &ione, x, &ione);      /* x := x + alpha*p         */
        F77_CALL(daxpy)(&n, &alpha, y, &ione, r, &ione);      /* r := r + alpha*A*p       */

        mfun(y, r, mdata);                      /* y := inv(M)r             */

        ry2 = ry;
        ry = F77_CALL(ddot)(&n, r, &ione, y, &ione);
        beta = ry / ry2;                        /* beta := r^Ty/r_^Ty_      */

        F77_CALL(dscal)(&n, &beta, p, &ione);
        F77_CALL(daxpy)(&n, &dminusone, y, &ione, p, &ione);       /* p = -y + beta*p          */

        //rnorm = F77_CALL(dnrm2(&n, r, &ione);
        //printf(" PCG> k = %4d,  r/b = %15.5e  (%15.5e ** %15.5e)\n", k, rnorm/bnorm, bnorm, tol);
    }
    free(r);
    free(y);
    free(p);

    /* pcg_status */
    if (k >= maxiter) {
        pcgstat->flag = PCG_MAXITER;
    } else {
        pcgstat->flag = PCG_OK;
    }
    pcgstat->relres   = rnorm/bnorm;
    pcgstat->iter     = k;
        
    return pcgstat->flag;
}

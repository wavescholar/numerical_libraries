/** \file   dmatrix.c
 *  \brief  Source file for matrix and vector manipulation functions.
 */

#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "def.h"
#include "dmatrix.h"

#ifndef Rpackage
    #include "blas.h"
    #include "lapack.h"
#else
    #include <R_ext/BLAS.h>
    #include <R_ext/Lapack.h>
#endif

#define PROFILING   FALSE

#if PROFILING
static long int time_yAx        = 0.0;
static long int time_yATx       = 0.0;
static long int time_yAmpqTx    = 0.0;
static long int time_yAmpqx     = 0.0;
static long int time_yHx        = 0.0;
static long int time_diagscale  = 0.0;
static long int time_baat       = 0.0;
static long int time_bata       = 0.0;
static long int time_aaxxtpa    = 0.0;
static long int time_potrs      = 0.0;
static long int time_posv       = 0.0;
#endif

/****************************************************************************/
/*                                                                          */
/*                           VECTOR OPERATIONS                              */
/*                                                                          */
/****************************************************************************/

/** \brief \f$ \|x\|_1 \f$
 *
 *  Returns 1-norm of a vector x.
 *
 *  @param  n   length of a vector x.
 *  @param  x   pointer to a vector x.
 *  @return     result.
 */
double dmat_norm1(const int n, const double *x)
{
    return F77_CALL(dasum)(&n, x, &ione);
}


/** \brief \f$ \|x\|_2 \f$
 *
 *  Returns 2-norm of a vector x.
 *
 *  @param  n   length of a vector x.
 *  @param  x   pointer to a vector x.
 *  @return     result.
 */
double dmat_norm2(const int n, const double *x)
{
    return F77_CALL(dnrm2)(&n, x, &ione);
}


/** \brief \f$ \|x\|_{\infty} \f$
 *
 *  Returns infinity-norm of a vector x
 *
 *  @param  n   length of a vector x
 *  @param  x   pointer to a vector x
 *  @return     result.
 */
double dmat_norminf(const int n, const double *x)
{
    return fabs(x[F77_CALL(idamax)(&n, x, &ione)-1]);
}


/** \brief \f$ x^Ty \f$
 *
 *  Returns dot product of a vector x and y.
 *
 *  @param  n   length of a vector x.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a vector y.
 *  @return     result.
 */
double dmat_dot(const int n, const double *x, const double *y)
{
    return F77_CALL(ddot)(&n, x, &ione, y, &ione);
}


/** \brief \f$ \mbox{dst}_i \leftarrow \mbox{val} \f$
 *
 *  Sets all the elements of a vector with a constant value.
 *
 *  @param  n   length of a vector.
 *  @param  val constant value to set.
 *  @param  dst pointer to a vector.
 */
void dmat_vset(int n, const double val, double *dst)
{
    while (n-- != 0)
        *dst++ = val;
}

void dmat_iset(int n, const int val, int *dst)
{
    while (n-- != 0)
        *dst++ = val;
}


/** \brief \f$ \mbox{dst} \leftarrow \mbox{src} \f$
 *
 *  Copies a vector.
 *
 *  @param  n   length of vectors.
 *  @param  src pointer to a source vector.
 *  @param  dst pointer to a destination vector.
 */
void dmat_vcopy(const int n, const double *src, double *dst)
{
    F77_CALL(dcopy)(&n, src, &ione, dst, &ione);
}

void dmat_icopy(const int n, const int *src, int *dst)
{
    memcpy(dst, src, sizeof(int)*n);
}


/** \brief \f$ y_i = \exp(x_i) \f$
 *
 *  Computes elementwise exp() of a vector.
 *
 *  @param  n   length of vectors.
 *  @param  x   pointer to a source vector.
 *  @param  y   pointer to a destination vector.
 */
void dmat_yexpx(const int n, const double *x, double *y)
{
#ifdef MKL 
    vdExp(n, x, y);
#else
    #ifdef ACML
        vrda_exp(n, x, y);
    #else
        const double *xi;
        double *yi;
        xi = x+n-1;
        yi = y+n-1;
        do {
            *yi-- = exp(*xi--);
        } while (xi >= x);
    #endif
#endif
}


/** \brief \f$ y_i = x_i^{1/2} \f$
 *
 *  Computes elementwise sqrt() of a vector.
 *
 *  @param  n   length of vectors.
 *  @param  x   pointer to a source vector.
 *  @param  y   pointer to a destination vector.
 */
void dmat_ysqrtx(const int n, const double *x, double *y)
{
#ifdef MKL 
    vdSqrt(n, x, y);
#else
    const double *xi;
    double *yi;
    xi = x+n-1;
    yi = y+n-1;
    do {
        *yi-- = sqrt(*xi--);
    } while (xi >= x);
#endif
}


/** \brief \f$ y_i = 1/x_i \f$
 *
 *  Computes elementwise inv() of a vector.
 *
 *  @param  n   length of vectors.
 *  @param  x   pointer to a source vector.
 *  @param  y   pointer to a destination vector.
 */
void dmat_yinvx(const int n, const double *x, double *y)
{
#ifdef MKL 
    vdInv(n, x, y);
#else
    const double *xi;
    double *yi;
    xi = x+n-1;
    yi = y+n-1;
    do {
        *yi-- = 1/(*xi--);
    } while (xi >= x);
#endif
}


/** \brief \f$ w = \alpha*x + \beta*y \f$
 *
 *  Computes weighted vector sum.
 *  - w = -x
 *  - w = alpha*x
 *  - w = -x + y
 *  - w = alpha*x + y
 *  - w = -x - y
 *  - w = alpha*x - y
 *  - w = alpha*x + beta*y
 *
 *  @param  n       length of vectors.
 *  @param  alpha   constant
 *  @param  x       pointer to a vector.
 *  @param  beta    constant
 *  @param  y       pointer to a vector.
 *  @param  w       pointer to a result vector.
 */
void dmat_waxpby(int n, double alpha, const double *x, double beta,
                 const double *y, double *w)

{
#if 1
    if (w != x && w != y)
    {
        dmat_vset(n, 0, w);
        F77_CALL(daxpy)(&n, &alpha, x, &ione, w, &ione);
        F77_CALL(daxpy)(&n, &beta , y, &ione, w, &ione);
    }
    else if (w == x && w == y)
    {
        double tmp;
        tmp = alpha+beta;
        F77_CALL(dscal)(&n, &tmp, w, &ione);
    }
    else if (w == x /*&& w != y */)
    {
        if (alpha != 1.0) F77_CALL(dscal)(&n, &alpha, w, &ione);
        if (beta  != 0.0) F77_CALL(daxpy)(&n, &beta , y, &ione, w, &ione);
    }
    else /* if (w == y && w != x ) */
    {
        if (beta  != 1.0) F77_CALL(dscal)(&n, &beta , w, &ione);
        if (alpha != 0.0) F77_CALL(daxpy)(&n, &alpha, x, &ione, w, &ione);
    }
#else
    int i;

    if (beta == 0.0)
    {
        if (alpha == -1.0)
        {
            for (i = 0; i < n; i++)
                w[i] = -x[i];
        }
        else
        {
            for (i = 0; i < n; i++)
                w[i] = alpha*x[i];
        }
    }
    else if (beta == 1.0)
    {
        if (alpha == -1.0)
        {
            for (i = 0; i < n; i++)
                w[i] = -x[i] + y[i];
        }
        else
        {
            for (i = 0; i < n; i++)
                w[i] = alpha*x[i] + y[i];
        }
    }
    else if (beta == -1.0)
    {
        if (alpha == -1.0)
        {
            for (i = 0; i < n; i++)
                w[i] = -x[i] - y[i];
        }
        else
        {
            for (i = 0; i < n; i++)
                w[i] = alpha*x[i] - y[i];
        }
    }
    else
    {
        for (i = 0; i < n; i++)
            w[i] = alpha*x[i] + beta*y[i];
    }
#endif
}


/**  \brief \f$ z_i = x_i*y_i \f$
 *
 *  Computes elementwise product of vectors.
 *  
 *  NOTE: x = x.*y is not allowed, i.e., w should not be the same with x.
 *
 *  @param  n   length of a vector x.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a vector y.
 *  @param  z   pointer to a result vector z.
 */
void dmat_elemprod(const int n, const double *x, const double *y, double *z)
{
    if (y != z) F77_CALL(dcopy)(&n, y, &ione, z, &ione);
    F77_CALL(dtbmv)("U", "N", "N", &n, &izero, x, &ione, z, &ione);
}

/** \brief \f$ z_i = x_i/y_i \f$
 *
 *  Computes elementwise division of vectors.
 *
 *  NOTE: y = x./y is not allowed, i.e., w should not be the same with y.
 *
 *  @param  n   length of a vector x.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a vector y.
 *  @param  z   pointer to a result vector z.
 */
void dmat_elemdivi(const int n, const double *x, const double *y, double *z)
{
#ifdef MKL
   vdDiv(n, x, y, z);
#else 
    if (x != z) F77_CALL(dcopy)(&n, x, &ione, z, &ione);
    F77_CALL(dtbsv)("U", "N", "N", &n, &izero, y, &ione, z, &ione);
#endif
}


/****************************************************************************/
/*                                                                          */
/*                           MATRIX OPERATIONS                              */
/*                                                                          */
/****************************************************************************/


/** \brief \f$ y = Ax \f$
 *
 *  Computes matrix-vector product.
 *
 *  @param  A   pointer to a matrix.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a result vector y.
 */
void dmat_yAx(const dmatrix *A, const double *x, double *y)
{
    int k, l;
    int m, n, nz;
    int *jdx, *rdx;
    double *val;
    double tmp;
#if PROFILING
    PROFILE_START
#endif

    m   = A->m;
    n   = A->n;
    nz  = A->nz;
    val = A->val;
    jdx = A->jdx;
    rdx = A->rdx;

    if (nz >= 0)    /* sparse matrix */
    {
        for (k = 0; k < m ; k++)
        {
            tmp = 0.0;
            for (l = rdx[k]; l < rdx[k+1]; l++)
            {
                tmp += val[l]*x[jdx[l]];
            }
            y[k] = tmp;
        }
    }
    else            /* dense matrix */
    {
        F77_CALL(dgemv)("T",&n,&m,&done,val,&n,x,&ione,&dzero,y,&ione);
    }

#if PROFILING
    PROFILE_END(time_yAx)
#endif
}


/** \brief \f$ y = A^Tx \f$
 *
 *  Computes transposed matrixr-vector product.
 *
 *  @param  A   pointer to a matrix.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a result vector y.
 */
void dmat_yATx(const dmatrix *A, const double *x, double *y)
{
    int k, l;
    int m, n, nz;
    int *jdx, *rdx;
    double *val;
    double tmp;
#if PROFILING
    PROFILE_START
#endif

    m   = A->m;
    n   = A->n;
    nz  = A->nz;
    val = A->val;
    jdx = A->jdx;
    rdx = A->rdx;

    dmat_vset(n, 0, y);
    if (nz >= 0)    /* sparse matrix */
    {
        for (k = 0; k < m ; k++)
        {
            tmp = x[k];
            for (l = rdx[k]; l < rdx[k+1]; l++)
            {
                y[jdx[l]] += val[l]*tmp;
            }
        }
    }
    else            /* dense matrix */
    {
        F77_CALL(dgemv)("N",&n,&m,&done,val,&n,x,&ione,&dzero,y,&ione);
    }
#if PROFILING
    PROFILE_END(time_yATx)
#endif
}

/** \brief \f$ y = (A-pq)x \f$
 *
 *  Computes row-rank updated matrix-vector product.
 *
 *  @param  A   pointer to a matrix.
 *  @param  p   pointer to a column vector.
 *  @param  q   pointer to a row vector.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a result vector y.
 */
void dmat_yAmpqx(const dmatrix *A, const double *p, const double *q,
                 const double *x, double *y)
{
    int k, l;
    int m, n, nz;
    int *jdx, *rdx;
    double *val;
    double tmp;
    double qx;
#if PROFILING
    PROFILE_START
#endif


    if (p == NULL) { dmat_yAx(A,x,y); return; }

    m   = A->m;
    n   = A->n;
    nz  = A->nz;
    val = A->val;
    jdx = A->jdx;
    rdx = A->rdx;

    qx = F77_CALL(ddot)(&n, q, &ione, x, &ione);

    if (nz >= 0)    /* sparse matrix */
    {
        for (k = 0; k < m ; k++)
        {
            tmp = 0.0;
            for (l = rdx[k]; l < rdx[k+1]; l++)
            {
                tmp += val[l]*x[jdx[l]];
            }
            y[k] = tmp - p[k]*qx;
        }
    }
    else
    {
        double mqx;
        mqx = -qx;
        F77_CALL(dcopy)(&m, p, &ione, y, &ione);
        F77_CALL(dgemv)("T",&n,&m,&done,val,&n,x,&ione,&mqx,y,&ione);
    }
#if PROFILING
    PROFILE_END(time_yAmpqx)
#endif
}


/** \brief \f$ y = (A-pq)^Tx \f$
 *
 *  Computes transposed row-rank updated matrix-vector product.
 *
 *  @param  A   pointer to a matrix.
 *  @param  p   pointer to a column vector.
 *  @param  q   pointer to a row vector.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a result vector y.
 */
void dmat_yAmpqTx(const dmatrix *A, const double *p, const double *q,
                 const double *x, double *y)
{
    int k, l;
    int m, n, nz;
    int *jdx, *rdx;
    double *val;
    double tmp;
    double px;
#if PROFILING
    PROFILE_START
#endif

    m   = A->m;
    n   = A->n;
    nz  = A->nz;
    val = A->val;
    jdx = A->jdx;
    rdx = A->rdx;

    if (p == NULL) {
        dmat_yATx(A,x,y);
        return;
    }

    dmat_vset(n, 0, y);
    px = -F77_CALL(ddot)(&m, p, &ione, x, &ione);
    F77_CALL(daxpy)(&n, &px, q, &ione, y, &ione);

    if (nz >= 0)    /* sparse matrix */
    {
        for (k = 0; k < m ; k++)
        {
            tmp = x[k];
            for (l = rdx[k]; l < rdx[k+1]; l++)
            {
                y[jdx[l]] += val[l]*tmp;
            }
        }
    }
    else
    {
        F77_CALL(dcopy)(&n, q, &ione, y, &ione);
        F77_CALL(dgemv)("N",&n,&m,&done,val,&n,x,&ione,&px,y,&ione);
    }

#if PROFILING
    PROFILE_END(time_yAmpqTx)
#endif
}


/** \brief \f$ y = A^TDAx \f$
 *
 *  Computes triple-matrix product.
 *
 *  @param  A   pointer to a matrix.
 *  @param  d   pointer to a diagonal vector \f$ D = diag(d) \f$
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a result vector y.
 */
void dmat_yATDAx(const dmatrix *A, const double *d, const double *x, double *y)
{
    int k, l;
    int m, n, nz;
    int *jdx, *rdx;
    double *val;
    double tmp;
    m   = A->m;
    n   = A->n;
    nz  = A->nz;
    val = A->val;
    jdx = A->jdx;
    rdx = A->rdx;

    if (nz >= 0)    /* sparse matrix */
    {
        for (k = 0; k < m ; k++)
        {
            tmp = 0.0;
            for (l = rdx[k]; l < rdx[k+1]; l++)
            {
                tmp += val[l]*x[jdx[l]];
            }
            tmp *= d[k];
            for (l = rdx[k]; l < rdx[k+1]; l++)
            {
                y[jdx[l]] = val[l]*tmp;
            }
        }
    }
    else
    {
        /* TO BE IMPLEMENTED */
        printf("dmat_yATDAx: dense case should be implemented\n");
        exit(0);
    }
}


/*           |b^T*D0*b  b^T*D0*A     0 | |v|
 *  y = Hx = |A^T*D0*b  A^T*D0*A+D1  D2|*|w|
 *           |0                  D2  D1| |u|
 */
/** \brief \f$
    y = Hx = \begin{array}{|ccc|}
            b^TD_0b & b^TD_0A    & 0   \\
            A^TD_0b & A^TD_0A+D1 & D_2 \\
            0       & D_2        & D_1 
            \end{array} \cdot \begin{array}{|c|}
            v \\ w \\ u
            \end{array}
\f$
 *
 *  Computes complex matrix-vector product.
 *
 *  @param  A   pointer to a matrix.
 *  @param  b   pointer to a vector.
 *  @param  d0  pointer to a vector.
 *  @param  d1  pointer to a vector.
 *  @param  d2  pointer to a vector.
 *  @param  x   pointer to a vector.
 *  @param  y   pointer to a result vector.
 */
void dmat_yHx_(const dmatrix *A, const double *b, const double *d0,
              const double *d1, const double *d2, const double *x, double *y)
{
    int k, l;
    int m, n, nz;
    int *jdx, *rdx;
    double *val;
    const double *v, *w, *u;
    double *vv, *ww, *uu;
    double tmp0, tmp1;

    m   = A->m;
    n   = A->n;
    nz  = A->nz;
    val = A->val;
    jdx = A->jdx;
    rdx = A->rdx;
    v   = &x[0]; w   = &x[1]; u   = &x[n+1];
    vv  = &y[0]; ww  = &y[1]; uu  = &y[n+1];

    if (nz >= 0)    /* sparse matrix */
    {
        /* tmp0 = b^T*D0*b*v */
        tmp0 = 0.0;
        for (k = 0; k < m ; k++)
        {
            tmp0 += b[k]*b[k]*d0[k];
        }
        tmp0 *= v[0];

        /* tmp0 += b^T*D0*A*w, ww = A^T*D0*[b A*w], uu = 0 */
        dmat_vset(n, 0, ww);
        dmat_vset(n, 0, uu);
        for (k = 0; k < m ; k++)
        {
            tmp1 = 0.0;
            for (l = rdx[k]; l < rdx[k+1]; l++)
            {
                tmp1 += val[l]*w[jdx[l]];
            }
            tmp0 += b[k]*d0[k]*tmp1;
            tmp1 = (tmp1+b[k]*v[0])*d0[k];

            for (l = rdx[k]; l < rdx[k+1]; l++)
            {
                ww[jdx[l]] += val[l]*tmp1;
            }
        }

        /* ww += D1*w+D2*u, uu += D2*w+D1*u */
        for (k = 0; k < n; k++)
        {
            ww[k] += d1[k]*w[k]+d2[k]*u[k];
            uu[k] += d2[k]*w[k]+d1[k]*u[k];
        }
        /* vv = tmp0 */
        vv[0] = tmp0;
    }
    else
    {
        /* TO BE IMPLEMENTED */
        printf("dmat_yATDAx: dense case should be implemented\n");
        exit(0);
    }
}

/** \brief \f$
    y = Hx = \begin{array}{|ccc|}
    b^TD_0b         & b^TD_0\tilde{A}            & 0   \\
    \tilde{A}^TD_0b & \tilde{A}^TD_0\tilde{A}+D1 & D_2 \\
    0               & D_2                        & D_1 
    \end{array} \cdot \begin{array}{|c|}
    v \\ w \\ u
    \end{array}
    \quad \mbox{  ,where   }
    \tilde{A} = A - pq
\f$
 *
 *  Computes complex matrix-vector product.
 *
 *  @param  A   pointer to a matrix.
 *  @param  p   pointer to a column vector.
 *  @param  q   pointer to a row vector.
 *  @param  b   pointer to a vector.
 *  @param  d0  pointer to a vector.
 *  @param  d1  pointer to a vector.
 *  @param  d2  pointer to a vector.
 *  @param  x   pointer to a vector.
 *  @param  y   pointer to a result vector.
 */
void dmat_yHx(const dmatrix *A, const double *p, const double *q,
              const double *b, const double *d0,
              const double *d1, const double *d2, const double *x, double *y)
{
    int k, l;
    int m, n, nz;
    int *jdx, *rdx;
    double *val;
    const double *v, *w, *u;
    double *vv, *ww, *uu;
    double tmp0, tmp1, tmp2;
    double qw;
#if PROFILING
    PROFILE_START
#endif


    if (p == NULL) 
    {
        dmat_yHx_(A,b,d0,d1,d2,x,y);
        return;
    }
    m   = A->m;
    n   = A->n;
    nz  = A->nz;
    val = A->val;
    jdx = A->jdx;
    rdx = A->rdx;
    v   = &x[0]; w   = &x[1]; u   = &x[n+1];
    vv  = &y[0]; ww  = &y[1]; uu  = &y[n+1];

    qw = F77_CALL(ddot)(&n, q, &ione, w, &ione);

    if (nz >= 0)    /* sparse matrix */
    {
        double mtmp2;
        /* tmp0 = b^T*D0*b*v */
        tmp0 = 0.0;
        tmp2 = 0.0;
        for (k = 0; k < m ; k++)
        {
            tmp0 += b[k]*b[k]*d0[k];
        }
        tmp0 *= v[0];

        /* tmp0 += b^T*D0*A*w, ww = A^T*D0*[b A*w] */
        dmat_vset(n, 0, ww);
        for (k = 0; k < m ; k++)
        {
            tmp1 = 0.0;
            for (l = rdx[k]; l < rdx[k+1]; l++)
            {
                tmp1 += val[l]*w[jdx[l]];
            }
            tmp1 -= p[k]*qw;

            tmp0 += b[k]*d0[k]*tmp1;
            tmp1 = (tmp1+b[k]*v[0])*d0[k];

            for (l = rdx[k]; l < rdx[k+1]; l++)
            {
                ww[jdx[l]] += val[l]*tmp1;
            }
            tmp2 += tmp1*p[k];
        }
        mtmp2 = -tmp2;
        F77_CALL(daxpy)(&n, &mtmp2, q, &ione, ww, &ione);

        /* ww += D1*w+D2*u, uu = D2*w+D1*u */
        for (k = 0; k < n; k++)
        {
            ww[k] += d1[k]*w[k]+d2[k]*u[k];
            uu[k] = d2[k]*w[k]+d1[k]*u[k];
        }
        /* vv = tmp0 */
        vv[0] = tmp0;
    }
    else
    {
        /* TO BE IMPLEMENTED */
        printf("dmat_yATDAx: dense case should be implemented\n");
        exit(0);
    }

#if PROFILING
    PROFILE_END(time_yHx)
#endif
}


/** \brief \f$ B_{ij} = A_{ij}^2 \f$
 *
 *  Computes elementwise matrix square.
 *
 *  @param  A   pointer to a matrix.
 *  @param  B   pointer to a result matrix.
 */
void dmat_elemAA(const dmatrix *A, dmatrix **B)
{
    dmatrix *tmp;
    int i, m, nz;
    double *val;
    m   = A->m;
    nz  = A->nz;
    val = A->val;

    tmp = malloc(sizeof(dmatrix));
    tmp->n  = A->n;
    tmp->m  = m;
    tmp->nz = nz;
    tmp->val = malloc(sizeof(double)*nz);
    tmp->idx = malloc(sizeof(int)*nz);
    tmp->jdx = malloc(sizeof(int)*nz);
    tmp->rdx = malloc(sizeof(int)*(m+1));
    memcpy(tmp->idx, A->idx, sizeof(int)*nz);
    memcpy(tmp->jdx, A->jdx, sizeof(int)*nz);
    memcpy(tmp->rdx, A->rdx, sizeof(int)*(m+1));
    for (i = 0; i < nz; i++)
    {
        tmp->val[i] = val[i]*val[i];
    }
    *B = tmp;
}


/** \brief \f$ M = D_L M D_R \f$
 *
 *  Computes left and right diagonal scaleing of M.
                    \f{eqnarray*}
                    M &=& D_L M \quad\mbox{ ,if dr = NULL }\\
                    M &=& M D_R \quad\mbox{ ,if dl = NULL }\\
                    M &=& D_L M D_R \quad\mbox{ ,otherwise }
                    \f}
 *
 *  @param  M       pointer to a matrix.
 *  @param  dl      pointer to a diagonal vector.
 *  @param  invl    inverse flag of left diagonal matrix.
                    \f{eqnarray*}
                    D_L &=& \mbox{diag}(dl)    \quad\mbox{ ,if invl} = 0\\
                    D_L &=& \mbox{diag}(1./dl) \quad\mbox{ ,if invl} = 1
                    \f}
 *  @param  dr      pointer to a diagonal vector.
 *  @param  invr    inverse flag of right diagonal matrix.
                    \f{eqnarray*}
                    D_R &=& \mbox{diag}(dr)    \quad\mbox{ ,if invr} = 0\\
                    D_R &=& \mbox{diag}(1./dr) \quad\mbox{ ,if invr} = 1
                    \f}
 */
void dmat_diagscale(dmatrix *M, const double *dl, const int invl,
                    const double *dr, const int invr)
{
    int i, j, m, n, nz;
    int *idx, *jdx;
    double *val;
#if PROFILING
    PROFILE_START
#endif

    m   = M->m;
    n   = M->n;
    nz  = M->nz;
    val = M->val;
    idx = M->idx;
    jdx = M->jdx;

    if (nz >= 0)    /* sparse matrix */
    {
        if (dl != NULL && dr != NULL)
        {
            if (invl == FALSE && invr == FALSE)
            {
                for (i = 0; i < nz; i++) val[i] *= dl[idx[i]]*dr[jdx[i]];
            }
            else if (invl == FALSE && invr == TRUE)
            {
                for (i = 0; i < nz; i++) val[i] *= dl[idx[i]]/dr[jdx[i]];
            }
            else if (invl == TRUE && invr == FALSE)
            {
                for (i = 0; i < nz; i++) val[i] *= dr[idx[i]]/dl[jdx[i]];
            }
            else
            {
                for (i = 0; i < nz; i++) val[i] /= dl[idx[i]]*dr[jdx[i]];
            }
        }
        else if (dl != NULL)
        {
            if (invl == FALSE)
            {
                for (i = 0; i < nz; i++) val[i] *= dl[idx[i]];
            }
            else
            {
                for (i = 0; i < nz; i++) val[i] /= dl[idx[i]];
            }
        }
        else if (dr != NULL)
        {
            if (invr == FALSE)
            {
                for (i = 0; i < nz; i++) val[i] *= dr[jdx[i]];
            }
            else
            {
                for (i = 0; i < nz; i++) val[i] /= dr[jdx[i]];
            }
        }
        else
        {
            /* no change */
        }
    }
    else    /* dense matrix */
    {
        if (dl != NULL && dr != NULL)
        {
            if (invl == FALSE && invr == FALSE)
            {
                for (i = 0; i < m; i++)
                    for (j = 0; j < n; j++)
                        val[i*n+j] *= dl[i]*dr[j];
            }
            else if (invl == FALSE && invr == TRUE)
            {
                for (i = 0; i < m; i++)
                    for (j = 0; j < n; j++)
                        val[i*n+j] *= dl[i]/dr[j];
            }
            else if (invl == TRUE && invr == FALSE)
            {
                for (i = 0; i < m; i++)
                    for (j = 0; j < n; j++)
                        val[i*n+j] *= dr[j]/dl[i];
            }
            else
            {
                for (i = 0; i < m; i++)
                    for (j = 0; j < n; j++)
                        val[i*n+j] /= dl[i]/dr[j];
            }
        }
        else if (dl != NULL)
        {
            if (invl == FALSE)
            {
                for (i = 0; i < m; i++)
                    for (j = 0; j < n; j++)
                        val[i*n+j] *= dl[i];
            }
            else
            {
                for (i = 0; i < m; i++)
                    for (j = 0; j < n; j++)
                        val[i*n+j] /= dl[i];
            }
        }
        else if (dr != NULL)
        {
            if (invr == FALSE)
            {
                for (i = 0; i < m; i++)
                    for (j = 0; j < n; j++)
                        val[i*n+j] *= dr[j];
            }
            else
            {
                for (i = 0; i < m; i++)
                    for (j = 0; j < n; j++)
                        val[i*n+j] /= dr[j];
            }
        }
        else
        {
            fprintf(stderr,"dmat_diagscale: wrong direction\n");
            exit(0);
        }
    }

#if PROFILING
    PROFILE_END(time_diagscale)
#endif
}


/** \brief \f$ M = M+\mbox{diag}(d) \f$
 *
 *  Adds diagonal entries to a matrix.
 *
 *  @param  M       pointer to a matrix.
 *  @param  d       pointer to a diagonal vector.
 */
void dmat_diagadd(dmatrix *M, const double *d)
{
    int m, n, nz;
    m  = M->m;
    n  = M->n;
    nz = M->nz;
    if (nz >= 0)
    {
        /* TO BE IMPLEMENTED */
        printf("dmat_diagadd: sparse case should be implemented\n");
        exit(0);
    }
    else
    {
        int i;
        double *val;
        val = M->val;
        for (i = 0; i < min(m,n); i++)
        {
            val[(n+1)*i] += d[i];
        }
    }
}


/** \brief \f$ y_j = \sum_{i} M_{ij} \f$
 *
 *  Computes column sum of a matrix.
 *
 *  @param  M       pointer to a matrix.
 *  @param  y       pointer to a result vector (column sum).
 */
void dmat_colsum(const dmatrix *M, double *y)
{
    int i, j, m, n, nz;
    int *idx, *jdx;
    double *val;
    n   = M->n;
    m   = M->m;
    nz  = M->nz;
    val = M->val;
    idx = M->idx;
    jdx = M->jdx;

    dmat_vset(n, 0, y);
    if (nz >= 0)
    {
        for (i = 0; i < nz; i++)
        {
            y[jdx[i]] += val[i];
        }
    }
    else    /* dense matrix */
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
            {
                y[j] += val[i*n+j];
            }
    }
}


/** \brief \f$ y_j = (1/m)\sum_{i} M_{ij} \f$
 *
 *  Computes column average of a matrix.
 *
 *  @param  M       pointer to a matrix.
 *  @param  y       pointer to a result vector (column average).
 */
void dmat_colavg(const dmatrix *M, double *y)
{
    int i, j, m, n, nz;
    int *idx, *jdx;
    double *val;
    n   = M->n;
    m   = M->m;
    nz  = M->nz;
    val = M->val;
    idx = M->idx;
    jdx = M->jdx;

    if (nz >= 0)
    {
        double tmp;
        dmat_vset(n, 0, y);
        for (i = 0; i < nz; i++)
        {
            y[jdx[i]] += val[i];
        }
        tmp = 1.0/m;
        F77_CALL(dscal)(&n, &tmp, y, &ione);
    }
    else    /* dense matrix */
    {
        for (j = 0; j < n; j++) /* col */
        {
            double sum = 0.0;
            for (i = 0; i < m; i++)
                sum += val[i*n+j];

            y[j] = sum/m;
        }
    }
}


/** \brief \f$ s_j = \left(\sum_{i}(M_{ij}-a_j)^2/(m-1)\right)^2 \f$
 *
 *  Computes column standard-deviation of a matrix.
 *
 *  @param  M       pointer to a matrix.
 *  @param  a       pointer to a column average vector.
 *  @param  s       pointer to a result vector (column std-dev).
 */
void dmat_colstd(const dmatrix *M, const double *a, double *s)
{
    int i, j, m, n, nz;
    int *jdx;
    double *val;
    m   = M->m;
    n   = M->n;
    nz  = M->nz;
    val = M->val;
    jdx = M->jdx;

    dmat_vset(n, 0, s);

    if (nz >= 0)
    {
        int *count;
        count = malloc(n*sizeof(int));
        for (i = 0; i < n; i++)
            count[i] = 0;

        for (i = 0; i < nz; i++)
        {
            j = jdx[i];
            s[j]  += (val[i]-a[j])*(val[i]-a[j]);
            count[j] += 1;
        }
        for (i = 0; i < n; i++)
        {
            s[i] = sqrt((s[i] + (m-count[i])*a[i]*a[i])/(m-1));
        }
        free(count);
    }
    else    /* dense matrix */
    {
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                s[j]  += (val[i*n+j]-a[j])*(val[i*n+j]-a[j]);
            }
        }
        for (i = 0; i < n; i++)
        {
            s[i] = sqrt(s[i]/(m-1));
        }
    }
}


/** \brief Duplicates a matrix.
 *
 *  Allocates all the memory and set sizes as in M, 
 *  but DO NOT copy the contents.
 *  To copy the contents, call dmat_copy.
 *
 *  @param  M       pointer to a source matrix.
 *  @param  dst     pointer to a destination matrix.
 */
void dmat_duplicate(const dmatrix* M, dmatrix** dst)
{
    int m, n, nz;
    dmatrix *mcp;

    mcp = malloc(sizeof(dmatrix));
    m   = M->m;
    n   = M->n;
    nz  = M->nz;
    mcp->m  = m;
    mcp->n  = n;
    mcp->nz = nz;
    if (nz >= 0)
    {
        mcp->val = malloc(sizeof(double)*nz);
        mcp->idx = malloc(sizeof(int)*nz);
        mcp->jdx = malloc(sizeof(int)*nz);
        mcp->rdx = malloc(sizeof(int)*(m+1));
    }
    else
    {
        mcp->val = malloc(sizeof(double)*(m*n));
        mcp->idx = NULL;
        mcp->jdx = NULL;
        mcp->rdx = NULL;
    }
    *dst = mcp;
}

/** \brief Copies the contents of a matrix.
 *
 *  Only copies the contents, NOT allocate memory.
 *
 *  @param  M       pointer to a source matrix.
 *  @param  dst     pointer to a destination matrix.
 */
void dmat_copy(const dmatrix* M, dmatrix* dst)
{
    int m, n, nz;

    m   = M->m;
    n   = M->n;
    nz  = M->nz;
    dst->m  = m;
    dst->n  = n;
    dst->nz = nz;
    if (nz >= 0)
    {
        memcpy(dst->val, M->val, sizeof(double)*nz);
        memcpy(dst->idx, M->idx, sizeof(int)*nz);
        memcpy(dst->jdx, M->jdx, sizeof(int)*nz);
        memcpy(dst->rdx, M->rdx, sizeof(int)*(m+1));
    }
    else
    {
        memcpy(dst->val, M->val, sizeof(double)*n*m);
        dst->idx = NULL;
        dst->jdx = NULL;
        dst->rdx = NULL;
    }
}


/** \brief Allocate memory for a dense matrix.
 *
 *  @param  M       pointer to an allocated matrix.
 *  @param  m       size of matrix (row).
 *  @param  n       size of matrix (column).
 */
void dmat_new_dense(dmatrix** M, const int m, const int n)
{
    dmatrix *dmat;
    dmat    = malloc(sizeof(dmatrix));
    dmat->m = m;
    dmat->n = n;
    dmat->nz = -1;
    dmat->val = malloc(sizeof(double)*m*n);
    dmat->idx = NULL;
    dmat->jdx = NULL;
    dmat->rdx = NULL;
    *M = dmat;
}

/** \brief Free memory for a matrix.
 *
 *  @param  M       pointer to a matrix to be freed.
 */
void dmat_free(dmatrix* M)
{
    if (M)
    {
        if (M->val) free(M->val);
        if (M->idx) free(M->idx);
        if (M->jdx) free(M->jdx);
        if (M->rdx) free(M->rdx);
        free(M);
    }
}


/** \brief Shows the content of a vector.
 *
 *  It only shows at most 10 elements of a vector.
 *
 *  @param  n       size of a vector.
 *  @param  v       pointer to a vector.
 */
void dmat_vprint(const int n, const double *v)
{
    int i;
    printf("\n");
    for (i = 0; i < min(n, 10); i++)
        printf("%5.4f    ",v[i]);
    printf("\n");
}


/** \brief Shows the content of a matrix.
 *
 *  It only shows at most 10 x 10 elements of a matrix.
 *
 *  @param  M       pointer to a matrix.
 */
void dmat_print(const dmatrix *M)
{
    int i, j, n, m, nz;
    int *idx, *jdx, *rdx;
    double *val;
    m   = M->m;
    n   = M->n;
    nz  = M->nz;
    val = M->val;
    idx = M->idx;
    jdx = M->jdx;
    rdx = M->rdx;

    if (nz >= 0)
    {
        double arr[10][10];
        for (i = 0; i < 10; i++)
            for (j = 0; j < 10; j++)
                arr[i][j] = 0;

        printf("\nsparse matrix: (%d x %d) nnz(%d)\n", m, n, nz);
        for (i = 0; i < nz; i++)
        {
            if (idx[i] < min(m,10) && jdx[i] < min(n,10))
                arr[idx[i]][jdx[i]] = val[i];
        }
        for (i = 0; i < min(m,10); i++)
        {
            for (j = 0; j < min(n,10); j++)
            {
                printf("%10g",arr[i][j]);
            }
            printf("\n");
        }
        /* print rdx */
        printf("-- rdx\n");
        for (i = 0; i < min(m+1,10); i++)
            printf("%d  ",rdx[i]);
        printf(" ... %d  %d\n", rdx[m-1], rdx[m]);
    }
    else
    {
        printf("\ndense matrix: (%d x %d) nnz(%d)\n", m, n, nz);
        for (i = 0; i < min(m,10); i++)
        {
            for (j = 0; j < min(n,10); j++)
            {
                printf("%10g",val[i*n+j]);
            }
            printf("\n");
        }
    }
}

void dmat_get_row(const dmatrix *M, const int rowidx, double *dst)
{
    if (M->nz >= 0)
    {
        int i;
        int *rdx, *jdx;
        double *val;

        rdx = M->rdx;
        jdx = M->jdx;
        val = M->val;

        dmat_vset(M->n, 0.0, dst);

        for (i = rdx[rowidx]; i < rdx[rowidx+1]; i++)
        {
            dst[jdx[i]] = val[i];
        }
    }
    else
    {
        F77_CALL(dcopy)(&(M->n), M->val+rowidx*M->n, &ione, dst, &ione);
    }
}


/** \brief Shows the summary of a matrix.
 *
 *  - sparsity.
 *  - size.
 *  - number of non-zero elements.
 *
 *  @param  M       pointer to a matrix.
 */
void dmat_summary(dmatrix *M)
{
    if (M->nz >= 0)
    {
        printf("Sparse matrix of size (%d x %d) with %d non-zeros\n",
                M->m, M->n, M->nz);
    }
    else
    {
        printf("Dense matrix of size (%d x %d)\n", M->m, M->n);
    }
}


/** \brief Build row index from csr info. (jdx and rdx)
 *
 *  - 
 *
 *  @param  M       pointer to a matrix.
 */
void dmat_build_idx(dmatrix *M)
{
    int i, j, m, nz;
    int *idx, *rdx;

    m = M->m;
    nz = M->nz;
    rdx = M->rdx;
    idx = M->idx;

    if (nz >= 0)
    {
        for (i = 0; i < m; i++)
        {
            for (j = rdx[i]; j < rdx[i+1]; j++) *idx++ = i;
        }
    }
    else
    {
        /* skip if dense matrix */
    }
}


/* BLAS LEVEL 2, 3 */

/** \brief \f$ B = AA^T \f$
 *
 *  Computes dense matrix-matrix-transpose product.
 *
 *  @param  A       pointer to a matrix.
 *  @param  B       pointer to a result matrix.
 */
void dmat_B_AAT(dmatrix *A, dmatrix *B)
{
#if PROFILING
    PROFILE_START
#endif

    F77_CALL(dsyrk)("L","T",&(A->m),&(A->n),&done,A->val,&(A->n),&dzero,B->val,&(B->n));

#if PROFILING
    PROFILE_END(time_baat)
#endif
}


/** \brief \f$ B = A^TA \f$
 *
 *  Computes dense matrix-transpose-matrix product.
 *
 *  @param  A       pointer to a matrix.
 *  @param  B       pointer to a result matrix.
 */
void dmat_B_ATA(dmatrix *A, dmatrix *B)
{
#if PROFILING
    PROFILE_START
#endif

    F77_CALL(dsyrk)("L","N",&(A->n),&(A->m),&done,A->val,&(A->n),&dzero,B->val,&(B->n));

#if PROFILING
    PROFILE_END(time_bata)
#endif
}


/** \brief \f$ A = \alpha x x^T +A \f$
 *
 *  Computes a row-rank update of a dense matrix.
 *
 *  @param  a       scalar.
 *  @param  x       pointer to a vector.
 *  @param  A       pointer to a matrix.
 */
void dmat_A_axxTpA(double a, double *x, dmatrix *A)
{
#if PROFILING
    PROFILE_START
#endif

    F77_CALL(dsyr)("L",&(A->n),&a,x,&ione,A->val,&(A->n));

#if PROFILING
    PROFILE_END(time_aaxxtpa)
#endif
}

/** \brief \f$ Ax = b \f$
 *
 *  Computes the solution of a linear system Ax = b via Cholesky method.
 *  Wrapper function to dports.
 *  That is, the cholesky factor is stored in A after computation.
 *
 *  @param  A       pointer to a matrix.
 *  @param  b       pointer to a vector.
 */
void dmat_potrs(const dmatrix *A, double *b)
{
#if PROFILING
    PROFILE_START
#endif

    int m = A->m;   /* A is symmetric, i.e.,  m = n */
    int one = 1;
    int info;
    F77_CALL(dpotrs)("L", &m, &one, A->val, &m, b, &m, &info);

#if PROFILING
    PROFILE_END(time_potrs)
#endif
}


/** \brief \f$ Ax = b \f$
 *
 *  Computes the solution of a linear system Ax = b via Cholesky method.
 *  Wrapper function to dposv.
 *  That is, A should be a Cholesky factor computed by potrs.
 *
 *  @param  A       pointer to a matrix (Cholesky factor).
 *  @param  b       pointer to a vector.
 */
void dmat_posv(const dmatrix *A, double *b)
{
#if PROFILING
    PROFILE_START
#endif

    int m = A->m;   /* A is symmetric, i.e.,  m = n */
    int one = 1;
    int info;
    F77_CALL(dposv)("L", &m, &one, A->val, &m, b, &m, &info);

#if PROFILING
    PROFILE_END(time_posv)
#endif
}


void dmat_profile()
{
#if PROFILING
printf(" time_yAx      %10ld\n", time_yAx);
printf(" time_yATx     %10ld\n", time_yATx);
printf(" time_yAmpqTx  %10ld\n", time_yAmpqTx);
printf(" time_yAmpqx   %10ld\n", time_yAmpqx);
printf(" time_yHx      %10ld\n", time_yHx);
printf(" time_diagscale%10ld\n", time_diagscale);
printf(" time_baat     %10ld\n", time_baat);
printf(" time_bata     %10ld\n", time_bata);
printf(" time_aaxxtpa  %10ld\n", time_aaxxtpa);
printf(" time_potrs    %10ld\n", time_potrs);
printf(" time_posv     %10ld\n", time_posv);
printf("\n total         %10ld\n", time_yAx+ time_yATx+ time_yAmpqTx+ time_yAmpqx+ time_yHx+ time_diagscale+ time_baat+ time_bata+ time_aaxxtpa+ time_potrs+ time_posv);
#endif
}



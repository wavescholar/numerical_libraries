/** \file   l1_logreg.c
 *  \brief  Main source file for l_1-regularized logistic regression
 *          problem solver.
 */

#if HAVE_CONFIG_H
#   include "config.h"
#endif

#define INTERNAL_PLOT 0

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "def.h"
#include "pcg.h"
#include "util.h"
#include "dmatrix.h"
#include "l1_logreg.h"

/* interior point method (IPM) parameters */
#define MU              2       /* updates t by a factor of MU          */
#define MAX_NT_ITER     400     /* terminates when iter > MAX_NT_ITER   */
#define ABSTOL          1e-8    /* terminates when duality gap < ABSTOL */

/* z optimizer parameters */
#define ZABSTOL         1e-12
#define MAX_NT_ITERZ    10
#define MAX_LS_ITERZ    10

/* line search parameters */
#define ALPHA           0.01
#define BETA            0.5
#define MAX_LS_ITER     100

/* pcg parameters */
#define MAX_PCG_ITER    5000

#define PRINT_BUF_SIZE  256

/*
 *  structure for problem data
 */
typedef struct
{
    dmatrix *matX1;
    dmatrix *matX2;
    double *ac;
    double *ar;

    double *b;
    double lambda;

    double *avg_x;
    double *std_x;
} problem_data_t;

void set_problem_data(problem_data_t * pdat, dmatrix * matX1, dmatrix * matX2,
                      double *ac, double *ar, double *b, double lambda,
                                                double *avg_x, double *std_x)
{
    pdat->matX1  = matX1;
    pdat->matX2  = matX2;
    pdat->ac     = ac;
    pdat->ar     = ar;
    pdat->b      = b;
    pdat->lambda = lambda;
    pdat->avg_x  = avg_x;
    pdat->std_x  = std_x;
}

void get_problem_data(problem_data_t * pdat, dmatrix ** matX1, dmatrix ** matX2,
                      double **ac, double **ar, double **b, double *lambda)
{
    *matX1  = pdat->matX1;
    *matX2  = pdat->matX2;
    *ac     = pdat->ac;
    *ar     = pdat->ar;
    *b      = pdat->b;
    *lambda = pdat->lambda;
}

void free_problem_data(problem_data_t *pdat)
{
    if (pdat->matX1) dmat_free(pdat->matX1);
    if (pdat->matX2) dmat_free(pdat->matX2);
    if (pdat->ac) free(pdat->ac);
    if (pdat->ar) free(pdat->ar);
    if (pdat->avg_x) free(pdat->avg_x);
    if (pdat->std_x) free(pdat->std_x);
    /* pdat->b shoud not be freed */
    /* pdat->lambda should not be freed */
}

/*
 *  structure for variables
 */
typedef struct
{
    double *x,  *v,  *w,  *u;
    double *dx, *dv, *dw, *du;
    double *gx, *gv, *gw, *gu;

    double *g;
    double *h;
    double *z;
    double *expz;
    double *expmz;

    double *d1;
    double *d2;
    double *Aw;

} variables_t;

void create_variables(variables_t * vars, int m, int n)
{
    int mbytes, nbytes, xbytes;

    nbytes = sizeof(double)*n;
    mbytes = sizeof(double)*m;
    xbytes = sizeof(double)*(1+n+n);

    /* allocate memory for vectors */

    vars->x  = malloc(xbytes);
    vars->v  = vars->x;
    vars->w  = vars->x+1;
    vars->u  = vars->x+1+n;

    vars->dx = malloc(xbytes);
    vars->dv = vars->dx;
    vars->dw = vars->dx+1;
    vars->du = vars->dx+1+n;

    vars->gx = malloc(xbytes);
    vars->gv = vars->gx;
    vars->gw = vars->gx+1;
    vars->gu = vars->gx+1+n;

    vars->g  = malloc(mbytes);
    vars->h  = malloc(mbytes);
    vars->z  = malloc(mbytes);
    vars->expz  = malloc(mbytes);
    vars->expmz = malloc(mbytes);

    vars->d1 = malloc(nbytes);
    vars->d2 = malloc(nbytes);
    vars->Aw = malloc(mbytes);

}

void free_variables(variables_t *vars)
{
    if (vars->x    ) free(vars->x);
    if (vars->dx   ) free(vars->dx);
    if (vars->gx   ) free(vars->gx);
    if (vars->g    ) free(vars->g);
    if (vars->h    ) free(vars->h);
    if (vars->z    ) free(vars->z);
    if (vars->expz ) free(vars->expz);
    if (vars->expmz) free(vars->expmz);
    if (vars->d1   ) free(vars->d1);
    if (vars->d2   ) free(vars->d2);
    if (vars->Aw   ) free(vars->Aw);
}

void get_variables(variables_t * vars, 
                   double **x,  double **v,  double **w,  double **u, 
                   double **dx, double **dv, double **dw, double **du,
                   double **gx, double **gv, double **gw, double **gu,
                   double **g,  double **h, 
                   double **z,  double **expz, double **expmz,
                   double **d1, double **d2,   double **Aw)
{
    *x  = vars->x;
    *v  = vars->v;
    *w  = vars->w;
    *u  = vars->u;
    *dx = vars->dx;
    *dv = vars->dv;
    *dw = vars->dw;
    *du = vars->du;
    *gx = vars->gx;
    *gv = vars->gv;
    *gw = vars->gw;
    *gu = vars->gu;
    *g  = vars->g;
    *h  = vars->h;
    *z  = vars->z;
    *expz  = vars->expz;
    *expmz = vars->expmz;
    *d1 = vars->d1;
    *d2 = vars->d2;
    *Aw = vars->Aw;
}

void allocate_temporaries(int m, int n, int is_pcg, double **tm1, double **tn1,
                          double **tn2, double **tn3, double **tn4,double **tx1,
                          double **precond, dmatrix **B, dmatrix **BB)
{
    int mbytes, nbytes, xbytes;

    nbytes = sizeof(double)*n;
    mbytes = sizeof(double)*m;
    xbytes = sizeof(double)*(n+n+1);

    *tm1 = malloc(mbytes);
    *tn1 = malloc(nbytes);

    /* problem specific storages */
    if (is_pcg)
    {
        *precond = malloc(sizeof(double)*(3*n+1));
        *B  = NULL;
        *BB = NULL;

        *tn2 = NULL;
        *tn3 = NULL;
        *tn4 = NULL;
        *tx1 = malloc(xbytes);
    }
    else
    {
        *precond = NULL;
        dmat_new_dense(B, m, n);
        dmat_new_dense(BB, min(m, n), min(m, n));
        if (n > m)
        {
            *tn2 = malloc(nbytes);
            *tn3 = malloc(nbytes);
            *tn4 = malloc(nbytes);
            *tx1 = malloc(xbytes);
        }
        else
        {
            *tn2 = malloc(nbytes);
            *tn3 = NULL;
            *tn4 = NULL;
            *tx1 = malloc(xbytes);
        }
    }
}

void free_temporaries(double *tm1, double *tn1, double *tn2, double *tn3,
                      double *tn4, double *tx1, double *precond,
                      dmatrix *B, dmatrix *BB)
{
    if (tm1) {free(tm1); tm1 = NULL;};
    if (tn1) {free(tn1); tn1 = NULL;};
    if (tn2) {free(tn2); tn2 = NULL;};
    if (tn3) {free(tn3); tn3 = NULL;};
    if (tn4) {free(tn4); tn4 = NULL;};
    if (tx1) {free(tx1); tx1 = NULL;};
    if (precond) free(precond);
    if (B ) dmat_free(B);
    if (BB) dmat_free(BB);
}


/*
 *  Data structures and functions for pcg
 */

/* structure for afun function in pcg */
typedef struct
{
    double *b;
    double *h;
    double *d1;
    double *d2;
    double *ac;
    double *ar;
    dmatrix *dmat;

} adata_t;

/* structure for mfun function in pcg */
typedef struct
{
    int m;
    int n;
    double *p0;
    double *p1;
    double *p2;
    double *p3;

} mdata_t;

void set_adata(adata_t * ad, dmatrix * dmat, double *ac, double *ar, 
               double *b, double *h, double *d1, double *d2)
{
    ad->dmat = dmat;
    ad->b    = b;
    ad->h    = h;
    ad->d1   = d1;
    ad->d2   = d2;
    ad->ac   = ac;
    ad->ar   = ar;
}

void set_mdata(mdata_t * md, int m, int n, double *precond)
{
    md->m  = m;
    md->n  = n;
    md->p0 = precond;
    md->p1 = precond+1;
    md->p2 = precond+1+n;
    md->p3 = precond+1+n+n;
}


/*
 *  Computes the logistic loss of a vector z, i.e.,
 *      sum(log(1+exp(-z))).
 *
 *  NOTE:
 *        To avoid numerical problems, we compute the loss as follows:
 *
 *          log(1+exp(-z))       =  log(1+exp(-z))      if (exp(z)<Inf)
 *                                  log(1+exp(+z))-z    otherwise
 */
# define LOG_INF    700                /* threshold value of x for exp(x)=INF */

/** \brief  Computes the logistic loss of vector z, i.e.,
 *          \f$ \sum_{i=1}^n \log(1 + \exp(-z_i)). \f$
 */
double logistic_loss(const int n, const double *z)
{
    int i;
    double ret = 0.0;

    for (i = 0; i < n; i++)
    {
        double mz = -z[i];
        ret += (mz < LOG_INF) ? log1p(exp(mz)) : log1p(exp(-mz)) + mz;
    }
    return ret;
}

/** \brief  Computes the logistic loss using precomputed exp(z), exp(-z), i.e.,
 *          \f$ \sum_{i=1}^n \log(1 + \exp(-z_i)). \f$
 */
double logistic_loss2(const int n, const double *z, 
                      const double *expz, const double *expmz)
{
    int i;
    double ret = 0.0;

    for (i = 0; i < n; i++)
    {
        ret += (z[i] > -LOG_INF) ? log1p(expmz[i]) : log1p(expz[i]) - z[i];
    }
    return ret;
}

/** \brief  Computes the sum of negative entropy of an n-vector z times n, i.e.,
 *          \f$ -\sum_{i=1}^n nz_i\log(nz_i)+(1-nz_i)log(1-nz_i). \f$
 */
double nentropy(const int n, const double *z)
{
    int i;
    double ret = 0.0;

    for (i = 0; i < n; i++)
    {
        double y1, y2;

        y1 = n * z[i];
        y2 = 1 - y1;

        /* handle 0*log(0) = 0 */
        ret -= (y1 > 0.0 && y2 > 0.0) ? y1 * log(y1) + y2 * log(y2) : 0.0;
    }
    return ret;
}

/** \brief  Computes the derivatives of the logistic loss.
 *          \f$ f'  = -1/(1+\exp(+z)) \f$,
 *          \f$ f'' = -1/(2+\exp(-z)+\exp(+z)) \f$.
 */
/*
 *      logistic_loss'   = -exp(z)/(1+exp(z))    = -1/(1+exp(+z))
 *      logistic_loss''  =  exp(z)/(1+exp(z))^2  =  1/((2+exp(z)+exp(-z))
 */
void fprimes(const int n, const double *expz, const double *expmz,
             double *f1prime, double *f2prime)
{
    int i;
    double ninv;

    ninv = 1.0 / n;

    for (i = 0; i < n; i++)
    {
        f1prime[i] = -ninv / (1.0 + expz[i]);
        f2prime[i] =  ninv / (2.0 + expz[i] + expmz[i]);
    }
}

/** \brief  Computes the derivatives of the logistic loss over v.
 */
void gradient_hessian_over_v(const int n, const double *expz, 
                             const double *expmz, const double *b, 
                             double *gradient, double *hessian)
{
    int i;
    double grad = 0.0;
    double hess = 0.0;

    for (i = 0; i < n; i++)
    {
        grad -= b[i] / (1.0 + expz[i]);
        hess += b[i] * b[i] / (2.0 + expmz[i] + expz[i]);
    }
    *gradient = grad;
    *hessian  = hess;
}

/** \brief  Evaluates the phi function.
 */
double eval_phi(const int m, const int n, const double *z, const double *expz,
                const double *expmz, const double *w, const double *u,
                const double lambda, const double t)
{
    int i;
    double penalty = 0.0;
    double barrier = 0.0;

    for (i = 0; i < n; i++)
    {
        penalty += u[i];
        barrier -= log( (-w[i] + u[i])*(w[i] + u[i]) );
    }
    return logistic_loss2(m,z,expz,expmz)/m + lambda*penalty + barrier/t;
}

/*
 *  Callback function of PCG, which computes Ax (in Ax = b).
 */
void afun(double *y, const double *x, const void *adata)
{
    adata_t *A;

    A = (adata_t *) adata;
    dmat_yHx(A->dmat, A->ac, A->ar, A->b, A->h, A->d1, A->d2, x, y);
}

/*
 *  Callback function of PCG, which computes M^{-1}y.
 */
void mfun(double *y, const double *x, const void *mdata)
{
    int i, n;
    const double *v, *w, *u;
    double *vv, *ww, *uu;
    double *p0, *p1, *p2, *p3;

    n  = ((mdata_t *) mdata)->n;
    p0 = ((mdata_t *) mdata)->p0;
    p1 = ((mdata_t *) mdata)->p1;
    p2 = ((mdata_t *) mdata)->p2;
    p3 = ((mdata_t *) mdata)->p3;

    v  = x;
    w  = x+1;
    u  = x+1+n;

    vv = y;
    ww = y+1;
    uu = y+1+n;

    vv[0] = v[0] / p0[0];
    for (i = 0; i < n; i++)
    {
        ww[i] =  p1[i]*w[i] - p2[i]*u[i];
        uu[i] = -p2[i]*w[i] + p3[i]*u[i];
    }
}

int is_indomain(double *w, double *u, int n)
{
    int i;
    int valid = 1;

    for (i = 0; i < n; i++)
    {
        valid &= (w[i] < u[i] && -w[i] < u[i]);
    }
    return valid;
}

/** \brief  Optimize intercept v.
 *
 */
void optimize_intercept(double *v, double *z, double *expz, double *expmz,
                        double *z2, const double *b, const double *Aw,
                        const int m)
{
    int i, j;

    for (i = 1; i <= MAX_NT_ITERZ; i++)
    {
        double g2, h2, nd2, dv2, s2, phi2;

        /* z = Aw + b*v */
        dmat_waxpby(m, v[0], b, 1, Aw, z);
        dmat_yexpx(m, z,    expz);
        dmat_yinvx(m, expz, expmz);

        gradient_hessian_over_v(m, expz, expmz, b, &g2, &h2);
        nd2 = g2 * g2 / h2 / 2;        /* newton decrement square */
        dv2 = -g2 / h2;                /* step */

        if (nd2 < ZABSTOL || i == MAX_NT_ITERZ)
        {
            break;
        }
        s2 = 1.0;
        phi2 = logistic_loss2(m, z, expz, expmz);

        for (j = 1; j < MAX_LS_ITERZ; j++) 
        {
            /* z2 = z + s2*dv2*b; */
            dmat_waxpby(m, s2 * dv2, b, 1, z, z2);

            if (logistic_loss(m, z2) - phi2 <= ALPHA * s2 * dv2 * g2) 
            {
                break;
            }
            s2 *= BETA;
        }
        v[0] += s2 * dv2;
    }
}

/** \brief  Compute search direction using pcg method.
 *
 */
void compute_searchdir_pcg(problem_data_t * pdat, variables_t * vars,
                           double t, double s, double gap, pcg_status_t * pcgstat,
                           adata_t * adata, mdata_t * mdata, double *precond,
                           double *tmp_m1, double *A2h, double *tmp_x1)
{
    int i, m, n, nz;
    double *p0, *p1, *p2, *p3;
    double normg, pcgtol, pcgmaxi, multfact;

    dmatrix *matX1, *matX2;
    double lambda, tinv;
    double *g, *h, *z, *expz, *expmz, *ac, *ar, *b, *d1, *d2, *Aw;
    double *x, *v, *w, *u, *dx, *dv, *dw, *du, *gv, *gw, *gu, *gx;

    static double pcgtol_factor = 1.0;


    get_problem_data(pdat, &matX1, &matX2, &ac, &ar, &b, &lambda);
    get_variables(vars, &x, &v, &w, &u, &dx, &dv, &dw, &du, &gx, &gv,
                  &gw, &gu, &g, &h, &z, &expz, &expmz, &d1, &d2, &Aw);
    m  = matX1->m;
    n  = matX1->n;
    nz = matX1->nz;
    tinv = 1.0 / t;

    p0 = &precond[0];
    p1 = &precond[1];
    p2 = &precond[1+n];
    p3 = &precond[1+n+n];

    /* dmat_vset(n+n+1, 0, dx); */

    dmat_yATx(matX2, h, A2h);        /* A2h = A2'*h */

    multfact = 0.0;
    if (ac != NULL)
    {
        /* h.*ac */
        dmat_elemprod(m, h, ac, tmp_m1);

        dmat_vset(n, 0, tmp_x1);
        dmat_yAmpqTx(matX1, NULL, NULL, tmp_m1, tmp_x1);
        dmat_elemprod(n, ar, tmp_x1, tmp_x1);

        for (i = 0; i < m; i++)
        {
            multfact += h[i] * ac[i] * ac[i];
        }
    }

    p0[0] = 0;
    for (i = 0; i < m; i++)
    {
        p0[0] += b[i] * b[i] * h[i];
    }

    /* complete forming gradient and d1, d2, precond */
    for (i = 0; i < n; i++)
    {
        double q1, q2, d3, div;

        q1 = 1.0 / (u[i] + w[i]);
        q2 = 1.0 / (u[i] - w[i]);

        gw[i] -= (q1 - q2) * tinv;        /* A'*g   - (q1-q2) */
        gu[i] = lambda - (q1 + q2) * tinv;        /* lambda - (q1+q2) */

        d1[i] = (q1 * q1 + q2 * q2) * tinv;
        d2[i] = (q1 * q1 - q2 * q2) * tinv;

        if (ac != NULL)
        {
            d3 = A2h[i] + d1[i] + multfact*ar[i]*ar[i] - 2*tmp_x1[i];
        }
        else
        {
            d3 = A2h[i] + d1[i];
        }
        div = 1 / (d3 * d1[i] - d2[i] * d2[i]);

        p1[i] = d1[i] * div;
        p2[i] = d2[i] * div;
        p3[i] = d3 * div;
    }
    normg = dmat_norm2(n+n+1, gx);

    pcgtol = min(1e-1, 0.3*gap/min(1.0,normg));
    /*
    pcgtol = min(1e-1, 0.3*gap/min(1.0,sqrt(normg)));
    */
    pcgmaxi = MAX_PCG_ITER;
    if (s < 1e-5)
    {
        pcgtol_factor *= 0.5;
    }
    else
    {
        pcgtol_factor = 1.0;
    }
     pcgtol = pcgtol*pcgtol_factor;

    dmat_waxpby(n+n+1, -1, gx, 0, NULL, tmp_x1);

    pcg(dx, pcgstat, afun, adata, mfun, mdata, tmp_x1, pcgtol, pcgmaxi, n+n+1);
}


/** \brief  Compute search direction using cholesky method (m < n).
 *
 */
void compute_searchdir_chol_fat(problem_data_t *pdat, variables_t *vars,
                                double t, dmatrix *B, dmatrix *BB,
                                double *tm1, double *bDA, 
                                double *d3inv, double *tmp31, double *tmp32)
{
    int i, m, n;
    double bDb;

    dmatrix *matX1, *matX2;
    double lambda, tinv;
    double *g, *h, *z, *expz, *expmz, *ac, *ar, *b, *d1, *d2, *Aw;
    double *x, *v, *w, *u, *dx, *dv, *dw, *du, *gv, *gw, *gu, *gx;

    get_problem_data(pdat, &matX1, &matX2, &ac, &ar, &b, &lambda);
    get_variables(vars, &x, &v, &w, &u, &dx, &dv, &dw, &du, &gx, &gv,
                  &gw, &gu, &g, &h, &z, &expz, &expmz, &d1, &d2, &Aw);
    m = matX1->m;
    n = matX1->n;
    tinv = 1.0 / t;

    /* bDb, Db */
    bDb = 0.0;
    for (i = 0; i < m; i++)
    {
        tm1[i] = h[i] * b[i];               /* tm1 = Db */
        bDb += b[i] * tm1[i];
    }

    /* bDA, D_inv */
    dmat_yATx(matX1, tm1, bDA);
    dmat_copy(matX1, B);                    /* B = A */
    dmat_yinvx(m, h, tm1);                  /* tm1 = D_inv */

    for (i = 0; i < n; i++)
    {
        double ui, wi, q1, q2, q3, gr2;

        ui = u[i];
        wi = w[i];
        q1 = 1.0 / (ui + wi);
        q2 = 1.0 / (ui - wi);
        q3 = ui * ui + wi * wi;

        gw[i] -= (q1 - q2) * tinv;          /* A'*g - (q1-q2)   */
        gu[i] = lambda - (q1 + q2)*tinv;    /* lambda - (q1+q2) */

        d1[i] = (q1 * q1 + q2 * q2) * tinv;
        d2[i] = (q1 * q1 - q2 * q2) * tinv;
        gr2 = gw[i] + 2 * gu[i] * ui * wi / q3;
        d3inv[i] = t * q3 / 2;              /* en = d3^{-1} */
        /* temporary use of tmp31 */
        tmp31[i] = sqrt(d3inv[i]);

        /* store temporary values in dw, du */
        dw[i] = d3inv[i] * bDA[i];          /* dw := d3inv.*bDA */
        du[i] = d3inv[i] * gr2;             /* du := d3inv.*gr2 */
    }
    /* B = B*D3^{1/2} */
    dmat_diagscale(B, NULL, FALSE, tmp31, FALSE);

    /* S = BB = ... */
    /* BB = A*D3_inv*A^T */
    dmat_B_AAT(B, BB);                      /* BB = B*B^T */

    /* BB = D_inv + A*D3_inv*A^T */
    dmat_diagadd(BB, tm1);

    /* SMW */
    dmat_yAx(matX1, dw, tm1);
    dmat_posv(BB, tm1);

    dmat_yATx(matX1, tm1, tmp31);
    dmat_elemprod(n, d3inv, tmp31, tmp31);
    dmat_waxpby(n, -1, tmp31, 1, dw, tmp31);

    dmat_yAx(matX1, du, tm1);

    dmat_potrs(BB, tm1);

    dmat_yATx(matX1, tm1, tmp32);
    dmat_elemprod(n, d3inv, tmp32, tmp32);
    dmat_waxpby(n, -1, tmp32, 1, du, tmp32);

    dv[0] = (-gv[0] + dmat_dot(n,bDA,tmp32)) / (bDb - dmat_dot(n,bDA,tmp31));

    /* dw = ... */
    dmat_waxpby(n, -dv[0], tmp31, -1, tmp32, dw);

    /* du = -(gu+d2.*dw)./d1; */
    for (i = 0; i < n; i++)
        du[i] = -(gu[i] + d2[i] * dw[i]) / d1[i];
}


/** \brief  Compute search direction using cholesky method (m > n).
 *
 */
void compute_searchdir_chol_thin(problem_data_t *pdat, variables_t *vars,
                                 double t, dmatrix *B, dmatrix *BB,
                                 double *tm1, double *bDA, double *d3)
{
    int i, m, n;
    double bDb, bDbinv;

    dmatrix *matX1, *matX2;
    double lambda, tinv;
    double *g, *h, *z, *expz, *expmz, *ac, *ar, *b, *d1, *d2, *Aw;
    double *x, *v, *w, *u, *dx, *dv, *dw, *du, *gv, *gw, *gu, *gx;

    get_problem_data(pdat, &matX1, &matX2, &ac, &ar, &b, &lambda);
    get_variables(vars, &x, &v, &w, &u, &dx, &dv, &dw, &du, &gx, &gv, 
                  &gw, &gu, &g, &h, &z, &expz, &expmz, &d1, &d2, &Aw);

    m = matX1->m;
    n = matX1->n;
    tinv = 1.0 / t;

    /* bDb, Db */
    bDb = 0.0;
    for (i = 0; i < m; i++)
    {
        tm1[i] = h[i] * b[i];        /* tm1 = Db */
        bDb += b[i] * tm1[i];
    }
    bDbinv = 1.0 / bDb;

    /* bDA */
    dmat_yATx(matX1, tm1, bDA);
    dmat_copy(matX1, B);                /* B = A */
    dmat_ysqrtx(m, h, tm1);        /* tm1 = D^{1/2} */

    /* B = D^{1/2}*B */
    dmat_diagscale(B, tm1, FALSE, NULL, FALSE);

    /* BB = A^T*D*A */
    dmat_B_ATA(B, BB);                /* BB = B^T*B */

    for (i = 0; i < n; i++)
    {
        double q1, q2, q3, ui, wi, gr2;

        ui = u[i];
        wi = w[i];
        q1 = 1.0 / (ui + wi);
        q2 = 1.0 / (ui - wi);
        q3 = ui * ui + wi * wi;

        gw[i] -= (q1 - q2) * tinv;        /* A'*g   - (q1-q2) */
        gu[i] = lambda - (q1 + q2) * tinv;        /* lambda - (q1+q2) */

        d1[i] = (q1 * q1 + q2 * q2) * tinv;
        d2[i] = (q1 * q1 - q2 * q2) * tinv;
        d3[i] = 2 / q3 * tinv;

        /*  dw = (bDA'*gv-bDb*gr2); */
        gr2 = gw[i] + 2 * gu[i] * ui * wi / q3;
        dw[i] = bDA[i] * gv[0] * bDbinv - gr2;
    }

    /* dw = (bDb*S-bDA'*bDA)\(bDA'*gv-bDb*gr2);
          = (S-1/bDb)*bDA'*bDA)\(bDA'*(gv/bDb)-gr2); */
    dmat_diagadd(BB, d3);
    dmat_A_axxTpA(-bDbinv, bDA, BB);
    dmat_posv(BB, dw);

    /* dv = (-bDA*dw-gv)/bDb; */
    dv[0] = -(dmat_dot(n, bDA, dw) + gv[0]) / bDb;

    /* du = -(gu+d2.*dw)./d1; */
    for (i = 0; i < n; i++)
        du[i] = -(gu[i] + d2[i] * dw[i]) / d1[i];
}


/** \brief  Perform backtraking linesearch.
 *
 */
double backtracking_linesearch(problem_data_t *pdat, variables_t *vars,
                               double t, double *Adw, double *xnew)
{
    int i, m, n;
    double phi, s;
    double *vnew, *wnew, *unew;

    dmatrix *matX1, *matX2;
    double *g, *h, *z, *expz, *expmz, *ac, *ar, *b, *d1, *d2, *Aw;
    double *x, *v, *w, *u, *dx, *dv, *dw, *du, *gv, *gw, *gu, *gx;
    double lambda, gdx;

    get_problem_data(pdat, &matX1, &matX2, &ac, &ar, &b, &lambda);
    m = matX1->m;
    n = matX1->n;
    get_variables(vars, &x, &v, &w, &u, &dx, &dv, &dw, &du, &gx, &gv,
                  &gw, &gu, &g, &h, &z, &expz, &expmz, &d1, &d2, &Aw);

    vnew = xnew;
    wnew = xnew+1;
    unew = xnew+1+n;

    s = 1.0;
    phi = eval_phi(m, n, z, expz, expmz, w, u, lambda, t);

    dmat_yAmpqx(matX1, ac, ar, dw, Adw);        /* see below */
    dmat_waxpby(m, dv[0], b, 1, Adw, Adw);      /* Adw := A*dw+b*dv */
    dmat_waxpby(m, v[0], b, 1, Aw, Aw);         /* Aw  := A*w+b*v   */
    dmat_waxpby(n+n+1, 1, dx, 1, x, xnew);      /* xnew:= x + dx    */
    gdx = dmat_dot(n+n+1, gx, dx);

    for (i = 1; i <= MAX_LS_ITER; i++)
    {
        /* x := x + s*dx */
        if (is_indomain(wnew, unew, n) == TRUE)
        {
            double newphi;

            /* z = A*(w+s*dw) + b*(v+s*dw)
                 = (Aw+bv) + s*(Adw+bdv),
               where Aw, Adw+bdv are vectors. */

            dmat_waxpby(m, s, Adw, 1, Aw, z);
            dmat_yexpx(m, z, expz);
            dmat_yinvx(m, expz, expmz);
            newphi = eval_phi(m, n, z, expz, expmz, wnew, unew, lambda, t);

            if (newphi <= phi + ALPHA * s * gdx) break;
        }
        s *= BETA;
        dmat_waxpby(n+n+1, s, dx, 1, x, xnew);
    }
    if (i > MAX_LS_ITER) return -1;

    dmat_vcopy(n+n+1,xnew,x);
    return s;
}

double backtracking_linesearch_deprecated(problem_data_t *pdat, 
                        variables_t *vars, double t, double *Adw, double *xnew)
{
    int i, m, n;
    double phi, s;

    dmatrix *matX1, *matX2;
    double *g, *h, *z, *expz, *expmz, *ac, *ar, *b, *d1, *d2, *Aw;
    double *x, *v, *w, *u, *dx, *dv, *dw, *du, *gv, *gw, *gu, *gx;
    double lambda, gdx;

    get_problem_data(pdat, &matX1, &matX2, &ac, &ar, &b, &lambda);
    m = matX1->m;
    n = matX1->n;
    get_variables(vars, &x, &v, &w, &u, &dx, &dv, &dw, &du, &gx, &gv,
                  &gw, &gu, &g, &h, &z, &expz, &expmz, &d1, &d2, &Aw);

    s = 1.0;
    phi = eval_phi(m, n, z, expz, expmz, w, u, lambda, t);

    dmat_yAmpqx(matX1, ac, ar, dw, Adw);        /* see below */
    dmat_waxpby(m, dv[0], b, 1, Adw, Adw);      /* Adw := A*dw+b*dv */
    dmat_waxpby(m, v[0], b, 1, Aw, Aw);         /* Aw  := A*w+b*v   */
    dmat_waxpby(n+n+1, 1, dx, 1, x, x);         /* x   := x + dx    */
    gdx = dmat_dot(n+n+1, gx, dx);

    for (i = 1; i <= MAX_LS_ITER; i++)
    {
        /* x := x + s*dx */

        if (is_indomain(w, u, n) == TRUE)
        {
            double newphi;

            /* z = A*(w+s*dw) + b*(v+s*dw)
                 = (Aw+bv) + s*(Adw+bdv),
               where Aw, Adw+bdv are vectors. */

            dmat_waxpby(m, s, Adw, 1, Aw, z);
            dmat_yexpx(m, z, expz);
            dmat_yinvx(m, expz, expmz);
            newphi = eval_phi(m, n, z, expz, expmz, w, u, lambda, t);

            if (newphi <= phi + ALPHA * s * gdx) break;
        }
        dmat_waxpby(n+n+1, -s*(1-BETA), dx, 1, x, x);
        s *= BETA;
    }
    if (i > MAX_LS_ITER) return -1;

    return s;
}

/*
 *  Prints progress message according to verbose level and mode (pcg/cholesky).
 */
void progress_pcg_v3(char *form, int in0, double in1, double in2, double in3, 
                     double in4, double in5, int in6, double in7, int in8) {
    fprintf(stderr, form, in0, in1, in2, in3, in4, in5, in6, in7, in8);
}
void progress_pcg_v2(char *form, int in0, double in1, double in2, double in3, 
                     double in4, double in5, int in6, double in7, int in8) {
    fprintf(stderr, form, in0, in1, in2, in8);
}
void progress_pcg_v0(char *form, int in0, double in1, double in2, double in3, 
                     double in4, double in5, int in6, double in7, int in8) {
    /* do nothing */
}
void progress_dir_v3(char *form, int in0, double in1, double in2, double in3, 
                     double in4, double in5, int in6, double in7, int in8) {
    fprintf(stderr, form, in0, in1, in2, in3, in4, in5);
}
void progress_dir_v2(char *form, int in0, double in1, double in2, double in3, 
                     double in4, double in5, int in6, double in7, int in8) {
    fprintf(stderr, form, in0, in1, in2);
}
void progress_dir_v0(char *form, int in0, double in1, double in2, double in3, 
                     double in4, double in5, int in6, double in7, int in8) {
    /* do nothing */
}

typedef void (*p2func_progress)(char *, int, double, double, double, double, double, int, double, int);

char *prt[] = {
    "    NT iter   ",
    "       gap    ",
    "    primal obj",
    "     dual obj ",
    "     step size",
    "     t value  ",
    "   pcg flag   ",
    "    pcg relres",
    "   pcg iter   ",
};
char *prn[] = {
    "%10d    ",
    "%14.4e",
    "%14.4e",
    "%14.4e",
    "%14.4e",
    "%14.4e",
    "%10d    ",
    "%14.4e",
    "%10d    ",
};
void init_progress(const int is_pcg, const int verbose_level,
                   char *format_buf, p2func_progress *print_progress)
{
    char name_buf[PRINT_BUF_SIZE];

    if (is_pcg) /* pcg */
    {
        switch (verbose_level)
        {
        case 3:
            sprintf(name_buf,"%s%s%s%s%s%s%s%s%s\n",prt[0],prt[1],
                prt[2],prt[3],prt[4],prt[5],prt[6],prt[7],prt[8]);
            sprintf(format_buf,"%s%s%s%s%s%s%s%s%s\n",prn[0],prn[1],
                prn[2],prn[3],prn[4],prn[5],prn[6],prn[7],prn[8]);
            *print_progress = &progress_pcg_v3;
            break;
        case 0:
        case 1:
            *print_progress = &progress_pcg_v0;
            break;
        default:
            /* default level 2 */
            sprintf(name_buf,"%s%s%s%s\n",prt[0],prt[1],prt[2],prt[8]);
            sprintf(format_buf,"%s%s%s%s\n",prn[0],prn[1],prn[2],prn[8]);
            *print_progress = &progress_pcg_v2;
        }
    }
    else /* direct */
    {
        switch (verbose_level)
        {
        case 3:
            sprintf(name_buf,"%s%s%s%s%s%s\n",prt[0],prt[1],
                prt[2],prt[3],prt[4],prt[5]);
            sprintf(format_buf,"%s%s%s%s%s%s\n",prn[0],prn[1],
                prn[2],prn[3],prn[4],prn[5]);
            *print_progress = &progress_dir_v3;
            break;
        case 0:
        case 1:
            *print_progress = &progress_dir_v0;
            break;
        default:
            /* default level 2 */
            sprintf(name_buf,"%s%s%s\n",prt[0],prt[1],prt[2]);
            sprintf(format_buf,"%s%s%s\n",prn[0],prn[1],prn[2]);
            *print_progress = &progress_dir_v2;
        }
    }
    if (verbose_level >= 2) printf(name_buf);
}

/******************************************************************************
 *                                                                            *
 *                            L1_LOGREG_TRAIN                                 *
 *                                                                            *
 ******************************************************************************/

/** \brief  Train (learn a model).
 *
 *  @param  X               train data (feature matrix)
 *  @param  b               train data (class vector)
 *  @param  lambda          regularization parameter
 *  @param  sflag           standardization flag
 *  @param  cflag           coefficients flag
 *  @param  verbose_level   verbose level
 *  @param  tolerance       duality gap (absolute) tolerance
 *  @param  initial_x       initial primal point (v0, w0, u0)
 *  @param  initial_t       initial barrier parameter
 *  @param  sol             model data (coefficients and intercept)
 *                          - sol[0]       : intercept
 *                          - sol[1..(n-1)] : coefficients
 *                          -- \f$ w./\sigma \f$ when standardized
 *                          -- \f$ w \f$ when not standardized
 *  @param  total_ntiter    number of newton iterations
 *  @param  total_pcgiter   number of pcg iterations
 */

int l1_logreg_train(dmatrix *X, double *b, double lambda, train_opts to,
                    double *initial_x, double *initial_t, double *sol,
                    int *total_ntiter, int *total_pcgiter)
{
    /* problem data */
    problem_data_t  prob;
    variables_t     vars;

    dmatrix *matX1;     /* matX1 = diag(b)*X_std */
    dmatrix *matX2;     /* matX2 = X_std.^2 (only for pcg) */
    double *ac, *ar;
    double *avg_x, *std_x;

    int m, n, ntiter, pcgiter, status;
    double pobj, dobj, gap;
    double t, s, maxAnu;

    double *g,  *h,  *z,  *expz, *expmz;
    double *x,  *v,  *w,  *u;
    double *dx, *dv, *dw, *du;
    double *gv, *gw, *gu, *gx;
    double *d1, *d2, *Aw;

    /* pcg variables */
    pcg_status_t pcgstat;
    adata_t adata;
    mdata_t mdata;
    double *precond;

    /* temporary variables */
    double *tm1, *tn1, *tn2, *tn3, *tn4, *tx1;

    /* temporary variables for dense case (cholesky) */
    dmatrix *B;     /* m x n (or m x n) */
    dmatrix *BB;    /* n x n (or m x m) */

    char format_buf[PRINT_BUF_SIZE];

#if INTERNAL_PLOT
    dmatrix *internal_plot;
    dmat_new_dense(&internal_plot, 3, MAX_NT_ITER);
    memset(internal_plot->val,0,sizeof(double)*3*MAX_NT_ITER);
    /* row 1: cum_nt_iter, row 2: cum_pcg_iter, row 3: duality gap */
#endif

    p2func_progress print_progress = NULL;

    /*
     *  INITIALIZATION
     */
    s       =  1.0;
    pobj    =  DBL_MAX;
    dobj    = -DBL_MAX;
    pcgiter =  0;
    matX1   =  NULL;
    matX2   =  NULL;

    init_pcg_status(&pcgstat);

    dmat_duplicate(X, &matX1);
    dmat_copy(X, matX1);

    m = matX1->m;
    n = matX1->n;

    if (to.sflag == TRUE)
    {
        /* standardize_data not only standardizes the data,
           but also multiplies diag(b). */
        standardize_data(matX1, b, &avg_x, &std_x, &ac, &ar);
    }
    else
    {
        /* matX1 = diag(b)*X */
        dmat_diagscale(matX1, b, FALSE, NULL, TRUE);
        avg_x = std_x = ac = ar = NULL;
    }

    if (matX1->nz >= 0)                /* only for pcg */
    {
        dmat_elemAA(matX1, &matX2);
    }
    else
    {
        matX2 = NULL;
    }

    set_problem_data(&prob, matX1, matX2, ac, ar, b, lambda, avg_x, std_x);

    create_variables(&vars, m, n);
    get_variables(&vars, &x, &v, &w, &u, &dx, &dv, &dw, &du, &gx, &gv,
                  &gw, &gu, &g, &h, &z, &expz, &expmz, &d1, &d2, &Aw);

    allocate_temporaries(m, n, (matX1->nz >= 0),
                         &tm1, &tn1, &tn2, &tn3, &tn4, &tx1, &precond, &B, &BB);

    if (initial_x == NULL)
    {
        dmat_vset(1, 0.0, v);
        dmat_vset(n, 0.0, w);
        dmat_vset(n, 1.0, u);
        dmat_vset(n+n+1, 0, dx);
        t = min(max(1.0, 1.0 / lambda), 2.0 * n / ABSTOL);
    }
    else
    {
        dmat_vcopy(n+n+1, initial_x, x);
        dmat_vset(n+n+1, 0, dx);
        t = *initial_t;
    }

    set_adata(&adata, matX1, ac, ar, b, h, d1, d2);
    set_mdata(&mdata, m, n, precond);

    /* select printing function and format according to
           verbose level and method type (pcg/direct) */

    if (to.verbose_level>=2) init_progress((matX1->nz >= 0), to.verbose_level,
                              format_buf, &print_progress);

    /*** MAIN LOOP ************************************************************/

    for (ntiter = 0; ntiter < MAX_NT_ITER; ntiter++)
    {
        /*
         *  Sets v as the optimal value of the intercept.
         */
        dmat_yAmpqx(matX1, ac, ar, w, Aw);
        optimize_intercept(v, z, expz, expmz, tm1, b, Aw, m);

        /*
         *  Constructs dual feasible point nu.
         */
        fprimes(m, expz, expmz, g, h);

        /* partially computes the gradient of phi.
           the rest part of the gradient will be completed while computing 
           the search direction. */

        gv[0] = dmat_dot(m, b, g);              /* gv = b'*g */
        dmat_yAmpqTx(matX1, ac, ar, g, gw);     /* gw = A'*g */

        dmat_waxpby(m, -1, g, 0, NULL, tm1);    /* nu = -g   */
        maxAnu = dmat_norminf(n, gw);           /* max(A'*nu) */

        if (maxAnu > lambda)
            dmat_waxpby(m, lambda / maxAnu, tm1, 0.0, NULL, tm1);

        /*
         *  Evaluates duality gap.
         */
        pobj = logistic_loss2(m,z,expz,expmz)/m + lambda*dmat_norm1(n,w);
        dobj = max(nentropy(m, tm1) / m, dobj);
        gap  = pobj - dobj;

#if INTERNAL_PLOT
        internal_plot->val[0*MAX_NT_ITER+ntiter] = (double)ntiter;
        internal_plot->val[1*MAX_NT_ITER+ntiter] = (double)pcgiter;
        internal_plot->val[2*MAX_NT_ITER+ntiter] = gap;
#endif
        if (to.verbose_level>=2)
        {
            (*print_progress)(format_buf, ntiter, gap, pobj, dobj, s, t,
                              pcgstat.flag, pcgstat.relres, pcgstat.iter);
        }

        /*
         *  Quits if gap < tolerance.
         */
        if (gap < to.tolerance ) /***********************************************/
        {
            if (sol != NULL)
            {
                /* trim solution */
                int i;
                double lambda_threshold;

                lambda_threshold = to.ktolerance*lambda;
                sol[0] = x[0];

                for (i = 0; i < n; i++)
                {
                    sol[i+1] = (fabs(gw[i])>lambda_threshold)? x[i+1] : 0.0;
                }
                /* if standardized, sol = coeff/std */
                if (to.sflag == TRUE && to.cflag == FALSE)
                {
                    dmat_elemdivi(n, sol+1, std_x, sol+1);
                    sol[0] -= dmat_dot(n, avg_x, sol+1);
                }
            }

            if (initial_x != NULL)
            {
                dmat_vcopy(n+n+1, x, initial_x);
                *initial_t = t;
            }

            if (total_pcgiter) *total_pcgiter = pcgiter;
            if (total_ntiter ) *total_ntiter  = ntiter;

            /* free memory */
            free_variables(&vars);
            free_temporaries(tm1, tn1, tn2, tn3, tn4, tx1, precond, B, BB);
            free_problem_data(&prob);

#if INTERNAL_PLOT
            write_mm_matrix("internal_plot",internal_plot,"",TYPE_G);
#endif
            return STATUS_SOLUTION_FOUND;

        } /********************************************************************/

        /*
         *  Updates t
         */
        if (s >= 0.5)
        {
            t = max(min(2.0 * n * MU / gap, MU * t), t);
        }
        else if (s < 1e-5)
        {
            t = 1.1*t;
        }

        /*
         *  Computes search direction.
         */
        if (matX1->nz >= 0)
        {
            /* pcg */
            compute_searchdir_pcg(&prob, &vars, t, s, gap, &pcgstat, &adata, 
                                  &mdata, precond, tm1, tn1, tx1);
            pcgiter += pcgstat.iter;
        }
        else
        {
            /* direct */
            if (n > m)
            {
                /* direct method for n > m, SMW */
                compute_searchdir_chol_fat(&prob, &vars, t, B, BB,
                                           tm1, tn1, tn2, tn3, tn4);
            }
            else
            {
                /* direct method for n <= m */
                compute_searchdir_chol_thin(&prob, &vars, t, B, BB,
                                            tm1, tn1, tn2);
            }
        }

        /*
         *  Backtracking linesearch & update x = (v,w,u) and z.
         */
        s = backtracking_linesearch(&prob, &vars, t, tm1, tx1);
        if (s < 0) break; /* BLS error */
    }
    /*** END OF MAIN LOOP *****************************************************/

    /*  Abnormal termination */
    if (s < 0)
    {
        status = STATUS_MAX_LS_ITER_EXCEEDED;
    }
    else /* if (ntiter == MAX_NT_ITER) */
    {
        status = STATUS_MAX_NT_ITER_EXCEEDED;
    }

    if (sol != NULL)
    {
        /* trim solution */
        int i;
        double lambda_threshold;

        lambda_threshold = to.ktolerance*lambda;
        sol[0] = x[0];

        for (i = 0; i < n; i++)
        {
            sol[i+1] = (fabs(gw[i])>lambda_threshold)? x[i+1] : 0.0;
        }
        /* if standardized, sol = coeff/std */
        if (to.sflag == TRUE && to.cflag == FALSE)
        {
            dmat_elemdivi(n, sol+1, std_x, sol+1);
            sol[0] -= dmat_dot(n, avg_x, sol+1);
        }
    }

    if (initial_x != NULL)
    {
        dmat_vcopy(n+n+1, x, initial_x);
        *initial_t = t;
    }
    if (total_pcgiter) *total_pcgiter = pcgiter;
    if (total_ntiter ) *total_ntiter  = ntiter;

    /* free memory */
    free_variables(&vars);
    free_temporaries(tm1, tn1, tn2, tn3, tn4, tx1, precond, B, BB);
    free_problem_data(&prob);

#if INTERNAL_PLOT
    write_mm_matrix("internal_plot",internal_plot,"",TYPE_G);
#endif
    return status;
}


/******************************************************************************
 *                                                                            *
 *                            L1_LOGREG_CLASSIFY                              *
 *                                                                            *
 ******************************************************************************/

/** \brief  Classify data.
 *
 *  \f$ \mbox{sign}\left(
        (X - 1\mu^T)\mbox{diag}(\sigma)^{-1}w +1v
    \right) =
        X\tilde{w} + 1(v-\mu^T\tilde{w})
    \f$
 *
 *  @param  X       test data (feature matrix)
 *  @param  b       test data (class vector)
 *  @param  sol     model data (coefficients and intercept)
 *                  - sol[0]       : intercept
 *                  - sol[1..(n-1)] : coefficients
                    -- \f$ w./\sigma \f$ when standardized
                    -- \f$ w \f$ when not standardized
 *  @param  result  test result vector
 *  @param  error_count numter of test error
 */
int l1_logreg_classify(dmatrix *X, double *b, double *sol,
                       int pflag, double *result, int *error_count)
{
    int m, n;

    m = X->m;
    n = X->n;

    if (result != NULL)
    {
        int     i;  /* example index */
        int     correct_count;

        correct_count = 0;

        dmat_yAx(X, sol+1, result);
        for (i = 0; i < m; i++)
            result[i] += sol[0];

        /* compare the estimation result 
            to the real class of examples if given */
        if (b != NULL)
        {
            for (i = 0; i < m; i++)
            {
                if ((result[i] > 0 && b[i] > 0) || 
                    (result[i] < 0 && b[i] < 0) ) {
                    result[i] = 1.0;    /* right */
                    correct_count++;
                } else {
                    result[i] = -1.0;   /* wrong */
                }
            }
            *error_count = m - correct_count;
        }
        else
        {
            /* classify test examples (works only when b==NULL) */
            if (pflag == TRUE)
            {
                for (i = 0; i < m; i++)
                    result[i] = 1.0/(1.0+exp(-result[i]));
            }
            else
            {
                for (i = 0; i < m; i++)
                    result[i] = (result[i]<0)? -1.0 : 1.0;
            }
                *error_count = -1;
        }
        return STATUS_OK;
    }
    else
    {
        *error_count = -1;
    }
    return STATUS_ERROR;
}






#ifndef PCG_H
#define PCG_H
/** \file   pcg.h
 *  \brief  Header file for a preconditioned conjugate gradient method
 *  implementation.
 */

#ifdef __cplusplus
extern "C" {
#endif

/** \struct pcg_status_t
 *  \brief  Contains pcg status information.
 */
typedef struct {
    int flag;           /**< exit flag; see below */
    int iter;           /**< number of iteration */
    double relres;      /**< relative residual */
} pcg_status_t;

/* defines for pcg_status.flag */
#define PCG_INIT        (-1)
#define PCG_OK          (0)
#define PCG_MAXITER     (1)

void init_pcg_status(pcg_status_t *pstst);

int pcg(double *x, pcg_status_t *pstat,
        void (*afun)(double*, const double*, const void*), const void* adata,
        void (*mfun)(double*, const double*, const void*), const void* mdata,
        double *b, double tol, double maxiter, int n);

#ifdef __cplusplus
}
#endif

#endif /* PCG_H */

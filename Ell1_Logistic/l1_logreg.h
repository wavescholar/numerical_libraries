#ifndef L1_LOGREG_H
#define L1_LOGREG_H

#include "dmatrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/* l1_logreg_train return values */
#define STATUS_SOLUTION_FOUND        0
#define STATUS_MAX_NT_ITER_EXCEEDED -1
#define STATUS_MAX_LS_ITER_EXCEEDED -2

typedef struct {
    int     sflag;          /* standardization flag */
    int     cflag;          /* c flag */
    int     verbose_level;  /* verbose level */
    double  tolerance;      /* tolerance */
    double  ktolerance;     /* kkt tolernace */
} train_opts;

int
l1_logreg_train(dmatrix *X, double *b, double lambda, train_opts to,
                double *initial_x, double *initial_t, double *sol,
                int *total_ntiter, int *total_pcgiter);

/* l1_logreg_classify return values */
#define STATUS_OK                   0
#define STATUS_ERROR               -1

int
l1_logreg_classify(dmatrix * matX, double *b, double *sol,
                   int pflag, double *result, int *error_count);


#ifdef __cplusplus
}
#endif

#endif /* L1_LOGREG_H */


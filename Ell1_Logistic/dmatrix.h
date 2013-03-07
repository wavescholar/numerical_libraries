#ifndef DMATRIX_H
#define DMATRIX_H
/** \file   dmatrix.h
 *  \brief  Header file for matrix and vector manipulation functions.
 */

#ifdef __cplusplus
extern "C" {
#endif

/** \struct dmatrix
 *  \brief  Structure for both dense and sparse matrix.
 *  - if the matrix is dense,
 *      then nz is set to -1 and idx, jdx are set to NULL.
 *  - if the matrix is sparse,
 *      then nz is set to the number of non-zero elements
 *      and the coordiante is saved at idx and jdx.
 *
 *  dmatrix stores a matrix in compressed sparse row (CSR) format.
 *  For computational efficiency, it also store the row indices
 *  as well as row start indices.
 *
 *  NOTE: As for sparse matrix, the matrix indices are not modified
 *  throughout the whole program. It is possible since
 *  matrix-matrix multiplication (except diagonal matrix) is
 *  not performed in PCG mode.
 */
typedef struct
{
    int     m;      /**< number of rows */
    int     n;      /**< number of columns */
    int     nz;     /**< number of non-zero entries*/
                    /**< nz >= 0: sparse, nz == -1: dense */
    double  *val;   /**< entry values */

    /* fields for sparse matrix */

    int     *jdx;   /**< column indices    (for both csr and coord) */
    int     *idx;   /**< row indices       (for coordinate) */
    int     *rdx;   /**< row start indices (for csr) */
    
} dmatrix;



double dmat_norm1(const int n, const double *x);
double dmat_norm2(const int n, const double *x);
double dmat_norminf(const int n, const double *x);

double dmat_dot(const int n, const double *x, const double *y);

void dmat_vset(int n, const double val, double *dst);
void dmat_vcopy(const int n, const double *src, double *dst);

void dmat_yexpx(const int n, const double *x, double *y);
void dmat_ysqrtx(const int n, const double *x, double *y);
void dmat_yinvx(const int n, const double *x, double *y);

void dmat_waxpby(int n, double alpha, const double *x, double beta,
                 const double *y, double *w);

void dmat_elemprod(const int n, const double *x, const double *y, double *z);
void dmat_elemdivi(const int n, const double *x, const double *y, double *z);



void dmat_yATx(const dmatrix *A, const double *x, double *y);

void dmat_yAx(const dmatrix *A, const double *x, double *y);

void dmat_yAmpqTx(const dmatrix *A, const double *p, const double *q,
                 const double *x, double *y);
void dmat_yAmpqx(const dmatrix *A, const double *p, const double *q,
                 const double *x, double *y);

void dmat_yHx_(const dmatrix *A, const double *b, const double *d0,
              const double *d1, const double *d2, const double *x, double *y);

void dmat_yHx(const dmatrix *A, const double *p, const double *q,
              const double *b, const double *d0,
              const double *d1, const double *d2, const double *x, double *y);

void dmat_elemAA(const dmatrix *A, dmatrix **A2);

void dmat_diagscale(dmatrix *M, const double *dl, const int invl,
                                const double *dr, const int invr);

void dmat_colsum(const dmatrix *M, double *vsum);

void dmat_colavg(const dmatrix *M, double *vavg);

void dmat_colstd(const dmatrix *M, const double *vavg, double *vstd);

void dmat_diagadd(dmatrix *M, const double *d);

void dmat_new_dense(dmatrix** M, const int m, const int n);

void dmat_free(dmatrix *M);

void dmat_duplicate(const dmatrix* M, dmatrix** Mcopy);
void dmat_copy(const dmatrix* M, dmatrix* dst);
void dmat_get_row(const dmatrix *M, const int rowidx, double *dst);

/* print for debugging */
void dmat_vprint(const int n, const double *v);

void dmat_print(const dmatrix *mat);
void dmat_summary(dmatrix *M);
void dmat_build_idx(dmatrix *M);

void dmat_B_AAT(dmatrix *A, dmatrix *B);
void dmat_B_ATA(dmatrix *A, dmatrix *B);
void dmat_A_axxTpA(double a, double *x, dmatrix *A);


void dmat_potrs(const dmatrix *A, double *b);
void dmat_posv(const dmatrix *A, double *b);

void dmat_profile();

#ifdef __cplusplus
}
#endif

#endif /* DMATRIX_H */

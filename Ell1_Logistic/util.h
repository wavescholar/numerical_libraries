#ifndef UTIL_H
#define UTIL_H
/** \file   util.h 
 *  \brief  Header file for utility functions.
 */

#include "dmatrix.h"

#ifdef __cplusplus
extern "C" {
#endif

#define     TYPE_E      0
#define     TYPE_G      1

/** \struct buffer
 *  \brief  simple buffer structure
 */
#define     BUFFER_INC_SIZE     1024
#define     BUFFER_SIZE_DOUBLE  sizeof(double)
#define     BUFFER_SIZE_INT     sizeof(int)
typedef struct {
    int elemsize;       /**< size of buffer element */
    int length;         /**< length of buffer */
    int curpos;         /**< size of buffer */
    void *ptr;          /**< puffer starting position */
} buffer_t;

void buffer_new(buffer_t **buf, int elemsize, int length);
void buffer_free(buffer_t *buf);
int  buffer_write(buffer_t *buf, void *src, int n);
int condense_solution(const int len, double *sol, int *idx);

/***/

double find_lambdamax(const dmatrix *X, const double *b, const int sflag);


void standardize_data(dmatrix *X, const double *b, double **average,
                      double **stddev, double **acol, double **arow);


int read_mm_new_matrix(const char *file_x, dmatrix **out_mat);
int read_mm_new_matrix_transpose(const char *file, dmatrix **out_mat);

int read_mm_new_vector(const char *file_x, double  **out_vec);

int write_mm_vector(const char *file, const int m, const double *vec,
                    const char *comments, const int type);
int write_mm_matrix(const char *file, dmatrix *mat,
                    const char *comments, const int type);

int write_mm_matrix_crd_header(FILE *fp, const int m, const int n,
                               const int nnz, const char *comments);
int write_mm_matrix_crd_column(FILE *fp, const int len,
                            const int col, const double *vec);
int write_mm_matrix_crd_column_threshold(FILE *fp, const int len,
                    const int col, const double *vec, const double threshold);
int get_mm_info(const char *file, int *m, int *n, int *nz);

int show_histogram(const int n, const double *coeff);

double userinput_threshold(void);
int thresholding(const int n, double *coeff, double threshold);


#ifdef __cplusplus
}
#endif

#endif /* UTIL_H */

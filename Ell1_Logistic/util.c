/** \file   util.c
 *  \brief  Source file for utility functions.
 */

#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>

#include "mmio.h"
#include "def.h"
#include "util.h"

#include "l1_logreg.h"
#include "dmatrix.h"


/** \brief Returns the maximum value of the regularization parameter lambda
 *         that gives a non-zero solution.
 *
 *  @param  X       feature matrix
 *  @param  b       class vector
 *  @param  sflag   standardization flag
 *                  - If sflag is 0, compute the maximum value of lambda
 *                  without standardization. 
 *                  - If sflag is 1, given matrix is standardized first and 
 *                  then the maximum value of lambda is computed.
 *
 *  @return         maximum value of lambda
 */
double find_lambdamax(const dmatrix *X, const double *b, const int sflag)
{
    double ret;
    int i, m, n;
    int mp, mn;
    double r1, r2;
    double *ar, *ac;
    double *tmp_m, *tmp_n;
    double *avg_x, *std_x;
    dmatrix *A;

    m = X->m;
    n = X->n;
    dmat_duplicate(X, &A);
    dmat_copy(X, A);

    tmp_m = malloc(m*sizeof(double));
    tmp_n = malloc(n*sizeof(double));

    if (sflag == TRUE)
    {
        standardize_data(A, b, &avg_x, &std_x, &ac, &ar);
    }
    else
    {
        dmat_diagscale(A, b, FALSE, NULL, TRUE);
        ac = ar = NULL;
    }

    /* number of positive class examples */
    mp = 0;
    for (i = 0; i < m; i++) 
    {
        mp += (b[i]>0 ? 1:0);
    }
    mn = m - mp;
    r1 = (double)mn/m;
    r2 = (double)mp/m;

    for (i = 0; i < m; i++)
    {
        tmp_m[i] = (b[i] > 0 ? r1 : r2);
    }
    dmat_yAmpqTx( A, ac, ar, tmp_m, tmp_n);
    ret = dmat_norminf(n, tmp_n) / m;

    free(tmp_n);
    free(tmp_m);
    dmat_free(A);
    if (avg_x) free(avg_x);
    if (std_x) free(std_x);
    if (ac) free(ac);
    if (ar) free(ar);

    return ret;
}


/** \brief Standardizes the data.
 *
 *  Standardizes the feature matrix.\n
 *  - If b is not given, compute \f$ X := (X-1^T\mu)\mbox{diag}(\sigma)^{-1} \f$
 *  - If b is given, compute
 *  \f$ X := \mbox{diag}(b)(X-1^T\mu)\mbox{diag}(\sigma)^{-1} \f$
 *  where \f$ \mu \f$ is column average and \f$ \sigma \f$ is column deviation.
 *
 *
 *  @param  X       feature matrix
 *  @param  b       class vector
 *  @param  average column average
 *  @param  stddev  column standard deviation
 *  @param  acol    column vector used for implicit standardization
 *  @param  arow    row vector used for implicit standardization
 */

void standardize_data(dmatrix *X, const double *b, double **average, 
                      double **stddev, double **acol, double **arow)
{
    int i, n, m, nz;
    double *avg, *std, *ac, *ar;

    n   = X->n;
    m   = X->m;
    nz  = X->nz;

    avg = malloc(n*sizeof(double));
    std = malloc(n*sizeof(double));
    dmat_colavg(X, avg);
    dmat_colstd(X, avg, std);

    for (i = 0; i < n; i++)
        if (std[i] < 1.0e-20) std[i] = 1;

    if ( nz >= 0)
    {
        ar  = malloc(n*sizeof(double));
        ac  = malloc(m*sizeof(double));

        /* X := diag(b)*X*inv(diag(std)) */
        dmat_diagscale(X, b, FALSE, std, TRUE);

        /* ac = diag(b)*1 */
        if (b != NULL) dmat_vcopy(m,   b, ac);
        else           dmat_vset (m, 1.0, ac);

        /* ar = avg^T*inv(diag(std)) */
        dmat_elemdivi(n, avg, std, ar);
    }
    else
    {
        ar  = NULL;
        ac  = NULL;

        for (i = 0; i < m; i++)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                X->val[i*n+j] -= avg[j];
            }
        }
        /* X = diag(b)*X*inv(diag(std)) */
        dmat_diagscale(X, b, FALSE, std, TRUE);
    }
    if (average != NULL) *average = avg; else free(avg);
    if (stddev  != NULL) *stddev  = std; else free(std);
    if (acol    != NULL) *acol    = ac;  else free(ac);
    if (arow    != NULL) *arow    = ar;  else free(ar);
}


/** \brief Read a matrix file.
 *
 *  Reads a Matrix Market formatted matrix from a file.
 *
 *  @param  file    matrix file name
 *  @param  out_mat matrix data
 *
 *  @return         result
 */

int read_mm_new_matrix(const char *file, dmatrix **out_mat)
{
    dmatrix *dmat;
    int i, j, m, n, nz;   
    FILE *fp;
    MM_typecode matcode;

    if ((fp = fopen(file, "r")) == NULL)
    {
        fprintf(stderr,"ERROR: Could not open file: %s\n",file);
        exit(1);
    }
    if (mm_read_banner(fp, &matcode) != 0) 
    {
        fprintf(stderr,"ERROR: Could not process Matrix Market banner.\n");
        exit(1);
    }
    if (!(mm_is_real(matcode) || mm_is_integer(matcode))) 
    {
        fprintf(stderr,"ERROR: Market Market type: [%s] not supported\n",
                mm_typecode_to_str(matcode));
        exit(1);
    }

    if (mm_is_sparse(matcode))
    {
        int cur_i, cur_j;
        double cur_val;
        double *val, *vtmp;
        int *idx, *jdx, *rdx;
        int *itmp, *jtmp, *tmp;

        if (mm_read_mtx_crd_size(fp, &m, &n, &nz) !=0) 
        {   /* find out size of sparse matrix */
            exit(1);
        }

        dmat = malloc(sizeof(dmatrix));

        itmp = malloc(nz*sizeof(int));
        jtmp = malloc(nz*sizeof(int));
        vtmp = malloc(nz*sizeof(double));

        rdx  = malloc((m+1)*sizeof(int));
        idx  = malloc(nz*sizeof(int));
        jdx  = malloc(nz*sizeof(int));
        val  = malloc(nz*sizeof(double));

        tmp  = malloc(m*sizeof(int));
        memset(tmp, 0, sizeof(int)*m);

        for (i = 0; i < nz; i++)
        {
            int ret;
            ret = fscanf(fp, "%d %d %lg\n", &cur_i, &cur_j, &cur_val);
            /* check format errors */
            if (ret != 3) {
                fprintf(stderr,"ERROR: Wrong file format, entry %d of %s.\n",
                        i+1, file);
                exit(1);
            }
            if (cur_i < 1 || cur_i > m || cur_j < 0 || cur_j > n ) {
                fprintf(stderr,"ERROR: Wrong index, entry %d of %s.\n",
                        i+1, file);
                exit(1);
            }
            itmp[i] = --cur_i;
            jtmp[i] = --cur_j;
            vtmp[i] = cur_val;
            tmp[ itmp[i] ]++ ;
        }

        rdx[0] = 0;
        for (i = 0; i < m ; i++)
        {
            rdx[i+1] = rdx[i] + tmp[i];
            tmp[i] = rdx[i];
        }
        for (i = 0; i < nz; i++)
        {
            int ii;
            ii = tmp[ itmp[i] ]++;
            idx[ii] = itmp[i];
            jdx[ii] = jtmp[i];
            val[ii] = vtmp[i];
        }
        dmat->n     = n;
        dmat->m     = m;
        dmat->nz    = nz;
        dmat->val   = val;
        dmat->idx   = idx;
        dmat->jdx   = jdx;
        dmat->rdx   = rdx;

        free(itmp);
        free(jtmp);
        free(vtmp);
        free(tmp);
    }
    else
    {
        double  *val;

        if (mm_read_mtx_array_size(fp, &m, &n) !=0) 
        {   /* find out size of dense matrix */
            exit(1);
        }

        dmat = malloc(sizeof(dmatrix));
        val  = malloc((m*n)*sizeof(double));

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                fscanf(fp, "%lg\n", &val[i*n+j]);
            }
        }
        dmat->n     = n;
        dmat->m     = m;
        dmat->nz    = -1;
        dmat->val   = val;
        dmat->idx   = NULL;
        dmat->jdx   = NULL;
        dmat->rdx   = NULL;
    }
    *out_mat = dmat;
    if (fp !=stdin) fclose(fp);

    return 0;
}

int read_mm_new_matrix_transpose(const char *file, dmatrix **out_mat)
{
    dmatrix *dmat;
    int i, j, m, n, nz;   
    FILE *fp;
    MM_typecode matcode;

    if ((fp = fopen(file, "r")) == NULL)
    {
        fprintf(stderr,"ERROR: Could not open file: %s\n",file);
        exit(1);
    }
    if (mm_read_banner(fp, &matcode) != 0) 
    {
        fprintf(stderr,"ERROR: Could not process Matrix Market banner.\n");
        exit(1);
    }
    if (!(mm_is_real(matcode) || mm_is_integer(matcode))) 
    {
        fprintf(stderr,"ERROR: Market Market type: [%s] not supported\n",
                mm_typecode_to_str(matcode));
        exit(1);
    }

    if (mm_is_sparse(matcode))
    {
        int cur_i, cur_j;
        double cur_val;
        double *val, *vtmp;
        int *idx, *jdx, *rdx;
        int *itmp, *jtmp, *tmp;

        if (mm_read_mtx_crd_size(fp, &n, &m, &nz) !=0) 
        {   /* find out size of sparse matrix */
            exit(1);
        }

        dmat = malloc(sizeof(dmatrix));

        itmp = malloc(nz*sizeof(int));
        jtmp = malloc(nz*sizeof(int));
        vtmp = malloc(nz*sizeof(double));

        rdx  = malloc((m+1)*sizeof(int));
        idx  = malloc(nz*sizeof(int));
        jdx  = malloc(nz*sizeof(int));
        val  = malloc(nz*sizeof(double));

        tmp  = malloc(m*sizeof(int));
        memset(tmp, 0, sizeof(int)*m);

        for (i = 0; i < nz; i++)
        {
            fscanf(fp, "%d %d %lg\n", &cur_j, &cur_i, &cur_val);
            itmp[i] = --cur_i;
            jtmp[i] = --cur_j;
            vtmp[i] = cur_val;
            tmp[ itmp[i] ]++ ;
        }

        rdx[0] = 0;
        for (i = 0; i < m ; i++)
        {
            rdx[i+1] = rdx[i] + tmp[i];
            tmp[i] = rdx[i];
        }
        for (i = 0; i < nz; i++)
        {
            int ii;
            ii = tmp[ itmp[i] ]++;
            idx[ii] = itmp[i];
            jdx[ii] = jtmp[i];
            val[ii] = vtmp[i];
        }
        dmat->n     = n;
        dmat->m     = m;
        dmat->nz    = nz;
        dmat->val   = val;
        dmat->idx   = idx;
        dmat->jdx   = jdx;
        dmat->rdx   = rdx;

        free(itmp);
        free(jtmp);
        free(vtmp);
        free(tmp);
    }
    else
    {
        double  *val;

        if (mm_read_mtx_array_size(fp, &n, &m) !=0) 
        {   /* find out size of dense matrix */
            exit(1);
        }

        dmat = malloc(sizeof(dmatrix));
        val  = malloc((m*n)*sizeof(double));

        for (j = 0; j < n*m; j++)
        {
            fscanf(fp, "%lg\n", &val[j]);
        }
        dmat->n     = n;
        dmat->m     = m;
        dmat->nz    = -1;
        dmat->val   = val;
        dmat->idx   = NULL;
        dmat->jdx   = NULL;
        dmat->rdx   = NULL;
    }
    *out_mat = dmat;
    if (fp !=stdin) fclose(fp);

    return 0;
}
/** \brief Read a vector file.
 *
 *  Reads a Matrix Market formatted vector from a file.
 *  Returns pointer to dense vector.
 *
 *  @param  file    vector file name
 *  @param  out_vec vector data
 *
 *  @return         result
 */
int read_mm_new_vector(const char *file, double **out_vec)
{
    double  *val;
    int i, m, n, nz;   
    FILE *fp;
    MM_typecode matcode;

    if ((fp = fopen(file, "r")) == NULL)
    {
        fprintf(stderr,"ERROR: Could not open file: %s\n",file);
        exit(1);
    }
    if (mm_read_banner(fp, &matcode) != 0) 
    {
        fprintf(stderr,"ERROR: Could not process Matrix Market banner.\n");
        exit(1);
    }
    if (!(mm_is_real(matcode) || mm_is_integer(matcode))) 
    {
        fprintf(stderr,"ERROR: Market Market type: [%s] not supported\n",
                mm_typecode_to_str(matcode));
        exit(1);
    }

    if (mm_is_sparse(matcode))
    {
        if (mm_read_mtx_crd_size(fp, &m, &n, &nz) !=0) 
        {   /* find out size of sparse matrix */
            exit(1);
        }
        if (n != 1)
        {
            fprintf(stderr,"ERROR: %s does not contain a column vector\n",file);
            exit(1);
        }

        val  = malloc(m*sizeof(double));
        memset(val, 0, m*sizeof(double));

        for (i = 0; i < nz; i++)
        {
            int cur_i, cur_j;
            double cur_val;

            fscanf(fp, "%d %d %lg\n", &cur_i, &cur_j, &cur_val);
            val[cur_i-1] = cur_val;
        }
    }
    else
    {
        if (mm_read_mtx_array_size(fp, &m, &n) !=0) 
        {   /* find out size of dense matrix */
            exit(1);
        }
        if (n != 1)
        {
            fprintf(stderr,"ERROR: %s does not contain a column vector\n",file);
            exit(1);
        }

        val  = malloc(m*sizeof(double));

        for (i = 0; i < m; i++)
        {
            double cur_val;
            fscanf(fp, "%lg\n", &cur_val);
            val[i] = cur_val;
        }
    }
    *out_vec = val;
    if (fp !=stdin) fclose(fp);

    return 0;
}


int write_mm_vector(const char *file, const int m, const double *vec,
                    const char *comments, const int type)
{
    int i;
    FILE *fp;
    MM_typecode matcode;

    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_array(&matcode);
    mm_set_real(&matcode);

    if ((fp = fopen(file, "w")) == NULL)
    {
        fprintf(stderr,"ERROR: Could not open file: %s\n",file);
        exit(1);
    }
    /* banner */
    mm_write_banner(fp, matcode); 
    /* comments */
    fprintf(fp, "%s", comments);
    /* size */
    mm_write_mtx_array_size(fp, m, 1);
    /* contents */

    if (type == TYPE_G)
    {
        for (i = 0; i < m; i++)
            fprintf(fp, "%10g\n", vec[i]);
    }
    else
    {
        for (i = 0; i < m; i++)
            fprintf(fp, "%22.15e\n", vec[i]);
    }


    fclose(fp);

    return 0;
}

int write_mm_matrix(const char *file, dmatrix *mat,
                    const char *comments, const int type)
{
	MM_typecode matcode;
    int i, j, m, n;
    FILE *fp;
    m = mat->m;
    n = mat->n;
    
    

    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);

    if (mat->nz < 0)
        mm_set_array(&matcode);
    else
        mm_set_coordinate(&matcode);

    mm_set_real(&matcode);

    if ((fp = fopen(file, "w")) == NULL)
    {
        fprintf(stderr,"ERROR: Could not open file: %s\n",file);
        exit(1);
    }
    /* banner */
    mm_write_banner(fp, matcode); 
    /* comments */
    fprintf(fp, "%s", comments);
    /* size */
    if (mat->nz < 0)
        mm_write_mtx_array_size(fp, m, n);
    else
        mm_write_mtx_crd_size(fp, m, n, mat->nz);

    /* contents */

    if (mat->nz < 0)
    {
        double *val;
        val = mat->val;

        if (type == TYPE_G)
        {
            for (j = 0; j < n; j++)
                for (i = 0; i < m; i++)
                    fprintf(fp, "%10g\n", val[i*n+j]);
        }
        else
        {
            for (j = 0; j < n; j++)
                for (i = 0; i < m; i++)
                    fprintf(fp, "%22.15e\n", val[i*n+j]);
        }
    }
    else
    {
        int *idx;
        int *jdx;
        double *val;
        idx = mat->idx;
        jdx = mat->jdx;
        val = mat->val;
        if (type == TYPE_G)
        {
            for (i = 0; i < mat->nz; i++)
            {
                fprintf(fp, "%d %d %15g\n", idx[i]+1, jdx[i]+1, val[i]);
            }
        }
        else
        {
            for (i = 0; i < mat->nz; i++)
            {
                fprintf(fp, "%d %d %22.15e\n", idx[i]+1, jdx[i]+1, val[i]);
            }
        }
    }
    fclose(fp);

    return 0;
}

int write_mm_matrix_crd_header(FILE *fp, const int m, const int n,
                               const int nnz, const char *comments)
{
    MM_typecode matcode;

    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);

    /* banner */
    mm_write_banner(fp, matcode); 
    /* comments */
    fprintf(fp, "%s", comments);
    /* size */
    mm_write_mtx_crd_size(fp, m, n, nnz);
    /* contents */

    return 0;
}

int write_mm_matrix_crd_column(FILE *fp, const int len,
                            const int col, const double *vec)
{
    int i;
    int count = 0;
    for (i = 0; i < len; i++)
        if ( vec[i] != 0 )
        {
            count++;
            fprintf(fp, "%d %d %22.15e\n", i+1, col+1, vec[i]);
        }

    return count;
}

int write_mm_matrix_crd_column_threshold(FILE *fp, const int len,
                    const int col, const double *vec, const double threshold)
{
    int i;
    int count = 0;
    for (i = 0; i < len; i++)
        if ( vec[i] < -threshold || threshold < vec[i] )
        {
            count++;
            fprintf(fp, "%d %d %22.15e\n", i+1, col+1, vec[i]);
        }

    return count;
}

int get_mm_info(const char *file, int *m, int *n, int *nz)
{
    FILE *fp;
    MM_typecode matcode;

    if ((fp = fopen(file, "r")) == NULL)
    {
        fprintf(stderr,"ERROR: Could not open file: %s\n",file);
        exit(1);
    }
    if (mm_read_banner(fp, &matcode) != 0) 
    {
        fprintf(stderr,"ERROR: Could not process Matrix Market banner.\n");
        exit(1);
    }
    if (!(mm_is_real(matcode) || mm_is_integer(matcode))) 
    {
        fprintf(stderr,"ERROR: Market Market type: [%s] not supported\n",
                mm_typecode_to_str(matcode));
        exit(1);
    }

    if (mm_is_sparse(matcode))
    {
        if (mm_read_mtx_crd_size(fp, m, n, nz) !=0) 
        {   /* find out size of sparse matrix */
            exit(1);
        }
    }
    else
    {
        if (mm_read_mtx_array_size(fp, m, n) !=0) 
        {   /* find out size of dense matrix */
            exit(1);
        }
        *nz = -1;
    }
    fclose(fp);

    return 0;
}

/** \brief  Shows histogram of coefficients in 10-base logscale.
 *
 *  @param  n       number of coefficients
 *  @param  coeff   vector of coefficients
 *
 */
#define BIN_EXP_BEGIN       (-7)
#define BIN_EXP_END         (+2)
#define BIN_NUM_PER_DECADE  (8)
#define BIN_IDX_NUM         ((BIN_EXP_END-BIN_EXP_BEGIN)*BIN_NUM_PER_DECADE)
#define BIN_IDX_BEGIN       (0)
#define BIN_IDX_END         (BIN_IDX_BEGIN+BIN_IDX_NUM-1)

int show_histogram(const int n, const double *coeff)
{
    int i, j;
    int bin[BIN_IDX_NUM];

    for (i = BIN_IDX_BEGIN; i <= BIN_IDX_END; i++)
        bin[i] = 0;

    for (i = 0; i < n; i++)
    {
        int val;
        val = (int)((log10(fabs(coeff[i]))-BIN_EXP_BEGIN)*BIN_NUM_PER_DECADE);
        if (val < BIN_IDX_BEGIN)
            val = BIN_IDX_BEGIN;
        else if (val > BIN_IDX_END)
            val = BIN_IDX_END;

        bin[val]++;
    }

    printf("Histogram of coefficients (10 base log scale, absolute values)\n");
    for (j = 5; j > 0; j--)
    {
        for (i = BIN_IDX_BEGIN; i <= BIN_IDX_END; i++)
            if (bin[i]>=j) printf("*");
            else           printf(" ");

        printf("\n");
    }
    for (i = BIN_EXP_BEGIN; i < BIN_EXP_END; i++)
    {
        int len;
        len = printf("|");
        for (j = 1; j <= BIN_NUM_PER_DECADE-len; j++)
            printf("-");
    }
    printf("\n");
    for (i = BIN_EXP_BEGIN; i < BIN_EXP_END; i++)
    {
        int len;
        len = printf("%d",i);
        for (j = 1; j <= BIN_NUM_PER_DECADE-len; j++)
            printf(" ");
    }
    printf("\n");

    return 0;
}

/** \brief  Gets a threshold value from standard input.
 *
 *  - If input value is positive, return it.
 *  - If input value is negative, assume it as a 10-based log value.
 *    For example, -3 is converted to 10^(-3).
 *
 *  @return             threshold value
 */
double userinput_threshold(void)
{
    float threshold;
    printf("\nInput threshold : ");
    scanf("%f",&threshold);
    if (threshold < 0)
        return pow(10, threshold);
    else
        return (double)threshold;
}

/** \brief  Thresholds coefficients.
 *
 *  Truncates coefficients if their absolute values are below the threshold.
 *
 *  @param n            number of coefficients
 *  @param coeff        vector of coefficients
 *  @param threshold    threshold value
 */
int thresholding(const int n, double *coeff, double threshold)
{
    int i;
    for (i = 0; i < n; i++)
    {
        if (-threshold< coeff[i] && coeff[i] < threshold)
            coeff[i] = 0.0; 
    }
    return 0;
}


/******************************/
void buffer_new(buffer_t **buf, int elemsize, int length)
{
    /* type is defined as the buffer element size */

    *buf = malloc(sizeof(buffer_t));
    (*buf)->ptr = malloc(length*elemsize);
    (*buf)->elemsize = elemsize;
    (*buf)->curpos = 0;
    (*buf)->length = length;
}

void buffer_free(buffer_t *buf)
{
    if (buf->ptr) free(buf->ptr);
    free(buf);
}

int buffer_write(buffer_t *buf, void *src, int n)
{
    void *ptr;
    ptr = buf->ptr;

    if (n > (buf->length - buf->curpos))
    {
        int alloc_len;
        alloc_len = (((buf->curpos+n)-1)/BUFFER_INC_SIZE+1)*BUFFER_INC_SIZE;

        ptr = realloc(ptr, alloc_len*buf->elemsize);
        if (ptr == NULL) return -1;
        buf->ptr = ptr;
        buf->length = alloc_len;
    }
	memcpy( (int*) buf->ptr , (int*)src, n*buf->elemsize );
    //memcpy( (buf->ptr+buf->curpos*buf->elemsize), (int*)src, n*(buf->elemsize );
	//memcpy( (buf->ptr+buf->curpos*buf->elemsize), src, n*buf->elemsize );
	//bbcrevisit changed this due to compile issue.  It looks like this method is only called 
	//when Rpackage is defined.

    buf->curpos += n;

    return 0;
}

/** \brief  Condense solution by removing zeros.
 *
 *  Truncates coefficients if their absolute values are below the threshold.
 *
 *  @param len          size of the original solution vector
 *  @param sol          solution vector
 *  @param idx          array of non-zero index in solution
 *
 *  ex: initial sol: [1 0 0 6 0 0 0 7 1]
 *      final sol:   [1 6 7 1]
 *            idx:   [1 4 8 9]
 */
int condense_solution(const int len, double *sol, int *idx)
{
    int i;
    int count = 0;

    int *idx_dst;
    double *sol_dst;
    idx_dst = idx;
    sol_dst = sol;

    for (i = 0; i < len; i++)
        if (sol[i] != 0)
        {
            count++;
            *sol_dst++ = sol[i];
            *idx_dst++ = i+1;   /* zero-base to one-base indexing */
        }
    
    return count;
}


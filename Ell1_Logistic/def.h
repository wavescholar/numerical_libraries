#ifndef DEF_H
#define DEF_H
/** \file   def.h
 *  \brief  Header file for global definitions
 */

#if HAVE_CONFIG_H
#   include "config.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* constant variables for Fortran blas arguments */
extern int     izero;
extern int     ione;
extern double  dzero;
extern double  done;
extern double  dminusone;

/* general definitions */
#define     TRUE                    (1)
#define     FALSE                   (0)


/* macro functions */
#define     max(x,y)                ((x)>(y)?(x):(y))
#define     min(x,y)                ((x)<(y)?(x):(y))


/* architecture dependent defs */
#if !HAVE_LOG1P
#   define log1p(x) log(1 + x)
#endif /* !HAVE_LOG1p */

#if !HAVE_GETTIMEOFDAY
#   define PROFILE_START
#   define PROFILE_END(x)

#else /* !HAVE_GETTIMEOFDAY */
#   define PROFILE_START                                \
    struct timeval start_time, end_time;                \
    gettimeofday(&start_time, NULL);
#   define PROFILE_END(x)                               \
    gettimeofday(&end_time, NULL);                      \
    x += (end_time.tv_sec-start_time.tv_sec) * 1000000 +\
                  (end_time.tv_usec-start_time.tv_usec);
#endif /* HAVE_GETTIMEOFDAY */



#ifdef __cplusplus
}
#endif

#endif /* DEF_H */

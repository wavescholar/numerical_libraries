
#ifndef _testablasunit_h
#define _testablasunit_h

#include "ap.h"
#include "ialglib.h"

#include "ablasf.h"
#include "ablas.h"


bool testablas(bool silent);


/*************************************************************************
Reference implementation

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
void refcmatrixrighttrsm(int m,
     int n,
     const ap::complex_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::complex_2d_array& x,
     int i2,
     int j2);


/*************************************************************************
Reference implementation

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
void refcmatrixlefttrsm(int m,
     int n,
     const ap::complex_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::complex_2d_array& x,
     int i2,
     int j2);


/*************************************************************************
Reference implementation

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
void refrmatrixrighttrsm(int m,
     int n,
     const ap::real_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::real_2d_array& x,
     int i2,
     int j2);


/*************************************************************************
Reference implementation

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
void refrmatrixlefttrsm(int m,
     int n,
     const ap::real_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::real_2d_array& x,
     int i2,
     int j2);


/*************************************************************************
Internal subroutine.
Triangular matrix inversion

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************/
bool internalcmatrixtrinverse(ap::complex_2d_array& a,
     int n,
     bool isupper,
     bool isunittriangular);


/*************************************************************************
Internal subroutine.
Triangular matrix inversion

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************/
bool internalrmatrixtrinverse(ap::real_2d_array& a,
     int n,
     bool isupper,
     bool isunittriangular);


/*************************************************************************
Reference SYRK subroutine.

  -- ALGLIB routine --
     16.12.2009
     Bochkanov Sergey
*************************************************************************/
void refcmatrixsyrk(int n,
     int k,
     double alpha,
     const ap::complex_2d_array& a,
     int ia,
     int ja,
     int optypea,
     double beta,
     ap::complex_2d_array& c,
     int ic,
     int jc,
     bool isupper);


/*************************************************************************
Reference SYRK subroutine.

  -- ALGLIB routine --
     16.12.2009
     Bochkanov Sergey
*************************************************************************/
void refrmatrixsyrk(int n,
     int k,
     double alpha,
     const ap::real_2d_array& a,
     int ia,
     int ja,
     int optypea,
     double beta,
     ap::real_2d_array& c,
     int ic,
     int jc,
     bool isupper);


/*************************************************************************
Reference GEMM,
ALGLIB subroutine
*************************************************************************/
void refcmatrixgemm(int m,
     int n,
     int k,
     ap::complex alpha,
     const ap::complex_2d_array& a,
     int ia,
     int ja,
     int optypea,
     const ap::complex_2d_array& b,
     int ib,
     int jb,
     int optypeb,
     ap::complex beta,
     ap::complex_2d_array& c,
     int ic,
     int jc);


/*************************************************************************
Reference GEMM,
ALGLIB subroutine
*************************************************************************/
void refrmatrixgemm(int m,
     int n,
     int k,
     double alpha,
     const ap::real_2d_array& a,
     int ia,
     int ja,
     int optypea,
     const ap::real_2d_array& b,
     int ib,
     int jb,
     int optypeb,
     double beta,
     ap::real_2d_array& c,
     int ic,
     int jc);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testablasunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testablasunit_test();


#endif


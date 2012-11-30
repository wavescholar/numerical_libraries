
#ifndef _testgq_h
#define _testgq_h

#include "ap.h"
#include "ialglib.h"

#include "hblas.h"
#include "reflections.h"
#include "creflections.h"
#include "sblas.h"
#include "ablasf.h"
#include "ablas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "hsschur.h"
#include "evd.h"
#include "gammafunc.h"
#include "gq.h"


/*************************************************************************
Test
*************************************************************************/
bool testgqunit(bool silent);


/*************************************************************************
Gauss-Hermite, another variant
*************************************************************************/
void buildgausshermitequadrature(int n,
     ap::real_1d_array& x,
     ap::real_1d_array& w);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testgq_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testgq_test();


#endif



#ifndef _testsvdunit_h
#define _testsvdunit_h

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
#include "bdsvd.h"
#include "svd.h"


/*************************************************************************
Testing SVD decomposition subroutine
*************************************************************************/
bool testsvd(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testsvdunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testsvdunit_test();


#endif


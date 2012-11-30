
#ifndef _testidwunit_h
#define _testidwunit_h

#include "ap.h"
#include "ialglib.h"

#include "tsort.h"
#include "nearestneighbor.h"
#include "reflections.h"
#include "hblas.h"
#include "creflections.h"
#include "sblas.h"
#include "ablasf.h"
#include "ablas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"
#include "hqrnd.h"
#include "matgen.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "xblas.h"
#include "densesolver.h"
#include "idwint.h"


/*************************************************************************
Testing IDW interpolation
*************************************************************************/
bool testidw(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testidwunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testidwunit_test();


#endif


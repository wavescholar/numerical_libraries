
#ifndef _testmlpeunit_h
#define _testmlpeunit_h

#include "ap.h"
#include "ialglib.h"

#include "mlpbase.h"
#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
#include "ablasf.h"
#include "ablas.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "matinv.h"
#include "linmin.h"
#include "minlbfgs.h"
#include "hblas.h"
#include "sblas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"
#include "xblas.h"
#include "densesolver.h"
#include "mlptrain.h"
#include "tsort.h"
#include "descriptivestatistics.h"
#include "bdss.h"
#include "mlpe.h"


bool testmlpe(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testmlpeunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testmlpeunit_test();


#endif


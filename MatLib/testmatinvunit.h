
#ifndef _testmatinvunit_h
#define _testmatinvunit_h

#include "ap.h"
#include "ialglib.h"

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


/*************************************************************************
Test
*************************************************************************/
bool testmatinv(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testmatinvunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testmatinvunit_test();


#endif


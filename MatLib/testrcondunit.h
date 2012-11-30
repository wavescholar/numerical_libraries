
#ifndef _testrcondunit_h
#define _testrcondunit_h

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


bool testrcond(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testrcondunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testrcondunit_test();


#endif


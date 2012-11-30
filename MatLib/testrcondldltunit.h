
#ifndef _testrcondldltunit_h
#define _testrcondldltunit_h

#include "ap.h"
#include "ialglib.h"

#include "ldlt.h"
#include "ssolve.h"
#include "estnorm.h"
#include "srcond.h"


bool testrcondldlt(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testrcondldltunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testrcondldltunit_test();


#endif


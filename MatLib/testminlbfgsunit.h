
#ifndef _testminlbfgsunit_h
#define _testminlbfgsunit_h

#include "ap.h"
#include "ialglib.h"

#include "linmin.h"
#include "minlbfgs.h"


bool testminlbfgs(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testminlbfgsunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testminlbfgsunit_test();


#endif


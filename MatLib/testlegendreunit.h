
#ifndef _testlegendreunit_h
#define _testlegendreunit_h

#include "ap.h"
#include "ialglib.h"

#include "legendre.h"


bool testlegendrecalculate(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testlegendreunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testlegendreunit_test();


#endif


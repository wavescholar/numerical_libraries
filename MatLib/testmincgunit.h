
#ifndef _testmincgunit_h
#define _testmincgunit_h

#include "ap.h"
#include "ialglib.h"

#include "linmin.h"
#include "mincg.h"


bool testmincg(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testmincgunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testmincgunit_test();


#endif


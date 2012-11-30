
#ifndef _testforestunit_h
#define _testforestunit_h

#include "ap.h"
#include "ialglib.h"

#include "tsort.h"
#include "descriptivestatistics.h"
#include "bdss.h"
#include "dforest.h"


bool testforest(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testforestunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testforestunit_test();


#endif


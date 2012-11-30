
#ifndef _testbdssunit_h
#define _testbdssunit_h

#include "ap.h"
#include "ialglib.h"

#include "tsort.h"
#include "descriptivestatistics.h"
#include "bdss.h"


/*************************************************************************
Testing BDSS operations
*************************************************************************/
bool testbdss(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testbdssunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testbdssunit_test();


#endif



#ifndef _testcorrunit_h
#define _testcorrunit_h

#include "ap.h"
#include "ialglib.h"

#include "ftbase.h"
#include "fft.h"
#include "conv.h"
#include "corr.h"


/*************************************************************************
Test
*************************************************************************/
bool testcorr(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testcorrunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testcorrunit_test();


#endif



#ifndef _testtrfacunit_h
#define _testtrfacunit_h

#include "ap.h"
#include "ialglib.h"

#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
#include "ablasf.h"
#include "ablas.h"
#include "trfac.h"


bool testtrfac(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testtrfacunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testtrfacunit_test();


#endif


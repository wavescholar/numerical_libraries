
#ifndef _testmatgenunit_h
#define _testmatgenunit_h

#include "ap.h"
#include "ialglib.h"

#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"


bool testmatgen(bool silent);


/*************************************************************************
Test whether matrix is SPD
*************************************************************************/
bool isspd(ap::real_2d_array a, int n, bool isupper);


bool obsoletesvddecomposition(ap::real_2d_array& a,
     int m,
     int n,
     ap::real_1d_array& w,
     ap::real_2d_array& v);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testmatgenunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testmatgenunit_test();


#endif


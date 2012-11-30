
#ifndef _testhqrndunit_h
#define _testhqrndunit_h

#include "ap.h"
#include "ialglib.h"

#include "hqrnd.h"


void calculatemv(const ap::real_1d_array& x,
     int n,
     double& mean,
     double& means,
     double& stddev,
     double& stddevs);


bool testhqrnd(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testhqrndunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testhqrndunit_test();


#endif


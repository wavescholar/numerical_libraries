
 
#include <stdio.h>
#include "testhermiteunit.h"

bool testhermitecalculate(bool silent)
{
    bool result;
    double err;
    double sumerr;
    double cerr;
    double threshold;
    int n;
    int maxn;
    int i;
    int j;
    int pass;
    ap::real_1d_array c;
    double x;
    double v;
    bool waserrors;

    err = 0;
    sumerr = 0;
    cerr = 0;
    threshold = 1.0E-9;
    waserrors = false;
    
    //
    // Testing Hermite polynomials
    //
    n = 0;
    err = ap::maxreal(err, fabs(hermitecalculate(n, double(1))-1));
    n = 1;
    err = ap::maxreal(err, fabs(hermitecalculate(n, double(1))-2));
    n = 2;
    err = ap::maxreal(err, fabs(hermitecalculate(n, double(1))-2));
    n = 3;
    err = ap::maxreal(err, fabs(hermitecalculate(n, double(1))+4));
    n = 4;
    err = ap::maxreal(err, fabs(hermitecalculate(n, double(1))+20));
    n = 5;
    err = ap::maxreal(err, fabs(hermitecalculate(n, double(1))+8));
    n = 6;
    err = ap::maxreal(err, fabs(hermitecalculate(n, double(1))-184));
    n = 7;
    err = ap::maxreal(err, fabs(hermitecalculate(n, double(1))-464));
    n = 11;
    err = ap::maxreal(err, fabs(hermitecalculate(n, double(1))-230848));
    n = 12;
    err = ap::maxreal(err, fabs(hermitecalculate(n, double(1))-280768));
    
    //
    // Testing Clenshaw summation
    //
    maxn = 10;
    c.setbounds(0, maxn);
    for(pass = 1; pass <= 10; pass++)
    {
        x = 2*ap::randomreal()-1;
        v = 0;
        for(n = 0; n <= maxn; n++)
        {
            c(n) = 2*ap::randomreal()-1;
            v = v+hermitecalculate(n, x)*c(n);
            sumerr = ap::maxreal(sumerr, fabs(v-hermitesum(c, n, x)));
        }
    }
    
    //
    // Testing coefficients
    //
    hermitecoefficients(0, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-1));
    hermitecoefficients(1, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-0));
    cerr = ap::maxreal(cerr, fabs(c(1)-2));
    hermitecoefficients(2, c);
    cerr = ap::maxreal(cerr, fabs(c(0)+2));
    cerr = ap::maxreal(cerr, fabs(c(1)-0));
    cerr = ap::maxreal(cerr, fabs(c(2)-4));
    hermitecoefficients(3, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-0));
    cerr = ap::maxreal(cerr, fabs(c(1)+12));
    cerr = ap::maxreal(cerr, fabs(c(2)-0));
    cerr = ap::maxreal(cerr, fabs(c(3)-8));
    hermitecoefficients(4, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-12));
    cerr = ap::maxreal(cerr, fabs(c(1)-0));
    cerr = ap::maxreal(cerr, fabs(c(2)+48));
    cerr = ap::maxreal(cerr, fabs(c(3)-0));
    cerr = ap::maxreal(cerr, fabs(c(4)-16));
    hermitecoefficients(5, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-0));
    cerr = ap::maxreal(cerr, fabs(c(1)-120));
    cerr = ap::maxreal(cerr, fabs(c(2)-0));
    cerr = ap::maxreal(cerr, fabs(c(3)+160));
    cerr = ap::maxreal(cerr, fabs(c(4)-0));
    cerr = ap::maxreal(cerr, fabs(c(5)-32));
    hermitecoefficients(6, c);
    cerr = ap::maxreal(cerr, fabs(c(0)+120));
    cerr = ap::maxreal(cerr, fabs(c(1)-0));
    cerr = ap::maxreal(cerr, fabs(c(2)-720));
    cerr = ap::maxreal(cerr, fabs(c(3)-0));
    cerr = ap::maxreal(cerr, fabs(c(4)+480));
    cerr = ap::maxreal(cerr, fabs(c(5)-0));
    cerr = ap::maxreal(cerr, fabs(c(6)-64));
    
    //
    // Reporting
    //
    waserrors = ap::fp_greater(err,threshold)||ap::fp_greater(sumerr,threshold)||ap::fp_greater(cerr,threshold);
    if( !silent )
    {
        printf("TESTING CALCULATION OF THE HERMITE POLYNOMIALS\n");
        printf("Max error                                 %5.2le\n",
            double(err));
        printf("Summation error                           %5.2le\n",
            double(sumerr));
        printf("Coefficients error                        %5.2le\n",
            double(cerr));
        printf("Threshold                                 %5.2le\n",
            double(threshold));
        if( !waserrors )
        {
            printf("TEST PASSED\n");
        }
        else
        {
            printf("TEST FAILED\n");
        }
    }
    result = !waserrors;
    return result;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testhermiteunit_test_silent()
{
    bool result;

    result = testhermitecalculate(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testhermiteunit_test()
{
    bool result;

    result = testhermitecalculate(false);
    return result;
}





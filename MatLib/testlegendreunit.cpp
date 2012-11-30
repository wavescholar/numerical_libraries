
 
#include <stdio.h>
#include "testlegendreunit.h"

bool testlegendrecalculate(bool silent)
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
    double t;
    bool waserrors;

    err = 0;
    sumerr = 0;
    cerr = 0;
    threshold = 1.0E-9;
    waserrors = false;
    
    //
    // Testing Legendre polynomials values
    //
    for(n = 0; n <= 10; n++)
    {
        legendrecoefficients(n, c);
        for(pass = 1; pass <= 10; pass++)
        {
            x = 2*ap::randomreal()-1;
            v = legendrecalculate(n, x);
            t = 1;
            for(i = 0; i <= n; i++)
            {
                v = v-c(i)*t;
                t = t*x;
            }
            err = ap::maxreal(err, fabs(v));
        }
    }
    
    //
    // Testing Clenshaw summation
    //
    maxn = 20;
    c.setbounds(0, maxn);
    for(pass = 1; pass <= 10; pass++)
    {
        x = 2*ap::randomreal()-1;
        v = 0;
        for(n = 0; n <= maxn; n++)
        {
            c(n) = 2*ap::randomreal()-1;
            v = v+legendrecalculate(n, x)*c(n);
            sumerr = ap::maxreal(sumerr, fabs(v-legendresum(c, n, x)));
        }
    }
    
    //
    // Testing coefficients
    //
    legendrecoefficients(0, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-1));
    legendrecoefficients(1, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-0));
    cerr = ap::maxreal(cerr, fabs(c(1)-1));
    legendrecoefficients(2, c);
    cerr = ap::maxreal(cerr, fabs(c(0)+double(1)/double(2)));
    cerr = ap::maxreal(cerr, fabs(c(1)-0));
    cerr = ap::maxreal(cerr, fabs(c(2)-double(3)/double(2)));
    legendrecoefficients(3, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-0));
    cerr = ap::maxreal(cerr, fabs(c(1)+double(3)/double(2)));
    cerr = ap::maxreal(cerr, fabs(c(2)-0));
    cerr = ap::maxreal(cerr, fabs(c(3)-double(5)/double(2)));
    legendrecoefficients(4, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-double(3)/double(8)));
    cerr = ap::maxreal(cerr, fabs(c(1)-0));
    cerr = ap::maxreal(cerr, fabs(c(2)+double(30)/double(8)));
    cerr = ap::maxreal(cerr, fabs(c(3)-0));
    cerr = ap::maxreal(cerr, fabs(c(4)-double(35)/double(8)));
    legendrecoefficients(9, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-0));
    cerr = ap::maxreal(cerr, fabs(c(1)-double(315)/double(128)));
    cerr = ap::maxreal(cerr, fabs(c(2)-0));
    cerr = ap::maxreal(cerr, fabs(c(3)+double(4620)/double(128)));
    cerr = ap::maxreal(cerr, fabs(c(4)-0));
    cerr = ap::maxreal(cerr, fabs(c(5)-double(18018)/double(128)));
    cerr = ap::maxreal(cerr, fabs(c(6)-0));
    cerr = ap::maxreal(cerr, fabs(c(7)+double(25740)/double(128)));
    cerr = ap::maxreal(cerr, fabs(c(8)-0));
    cerr = ap::maxreal(cerr, fabs(c(9)-double(12155)/double(128)));
    
    //
    // Reporting
    //
    waserrors = ap::fp_greater(err,threshold)||ap::fp_greater(sumerr,threshold)||ap::fp_greater(cerr,threshold);
    if( !silent )
    {
        printf("TESTING CALCULATION OF THE LEGENDRE POLYNOMIALS\n");
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
bool testlegendreunit_test_silent()
{
    bool result;

    result = testlegendrecalculate(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testlegendreunit_test()
{
    bool result;

    result = testlegendrecalculate(false);
    return result;
}






 
#include <stdio.h>
#include "testlaguerreunit.h"

bool testlaguerrecalculate(bool silent)
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
    // Testing Laguerre polynomials
    //
    n = 0;
    err = ap::maxreal(err, fabs(laguerrecalculate(n, 0.5)-1.0000000000));
    n = 1;
    err = ap::maxreal(err, fabs(laguerrecalculate(n, 0.5)-0.5000000000));
    n = 2;
    err = ap::maxreal(err, fabs(laguerrecalculate(n, 0.5)-0.1250000000));
    n = 3;
    err = ap::maxreal(err, fabs(laguerrecalculate(n, 0.5)+0.1458333333));
    n = 4;
    err = ap::maxreal(err, fabs(laguerrecalculate(n, 0.5)+0.3307291667));
    n = 5;
    err = ap::maxreal(err, fabs(laguerrecalculate(n, 0.5)+0.4455729167));
    n = 6;
    err = ap::maxreal(err, fabs(laguerrecalculate(n, 0.5)+0.5041449653));
    n = 7;
    err = ap::maxreal(err, fabs(laguerrecalculate(n, 0.5)+0.5183392237));
    n = 8;
    err = ap::maxreal(err, fabs(laguerrecalculate(n, 0.5)+0.4983629984));
    n = 9;
    err = ap::maxreal(err, fabs(laguerrecalculate(n, 0.5)+0.4529195204));
    n = 10;
    err = ap::maxreal(err, fabs(laguerrecalculate(n, 0.5)+0.3893744141));
    n = 11;
    err = ap::maxreal(err, fabs(laguerrecalculate(n, 0.5)+0.3139072988));
    n = 12;
    err = ap::maxreal(err, fabs(laguerrecalculate(n, 0.5)+0.2316496389));
    
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
            v = v+laguerrecalculate(n, x)*c(n);
            sumerr = ap::maxreal(sumerr, fabs(v-laguerresum(c, n, x)));
        }
    }
    
    //
    // Testing coefficients
    //
    laguerrecoefficients(0, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-1));
    laguerrecoefficients(1, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-1));
    cerr = ap::maxreal(cerr, fabs(c(1)+1));
    laguerrecoefficients(2, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-double(2)/double(2)));
    cerr = ap::maxreal(cerr, fabs(c(1)+double(4)/double(2)));
    cerr = ap::maxreal(cerr, fabs(c(2)-double(1)/double(2)));
    laguerrecoefficients(3, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-double(6)/double(6)));
    cerr = ap::maxreal(cerr, fabs(c(1)+double(18)/double(6)));
    cerr = ap::maxreal(cerr, fabs(c(2)-double(9)/double(6)));
    cerr = ap::maxreal(cerr, fabs(c(3)+double(1)/double(6)));
    laguerrecoefficients(4, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-double(24)/double(24)));
    cerr = ap::maxreal(cerr, fabs(c(1)+double(96)/double(24)));
    cerr = ap::maxreal(cerr, fabs(c(2)-double(72)/double(24)));
    cerr = ap::maxreal(cerr, fabs(c(3)+double(16)/double(24)));
    cerr = ap::maxreal(cerr, fabs(c(4)-double(1)/double(24)));
    laguerrecoefficients(5, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-double(120)/double(120)));
    cerr = ap::maxreal(cerr, fabs(c(1)+double(600)/double(120)));
    cerr = ap::maxreal(cerr, fabs(c(2)-double(600)/double(120)));
    cerr = ap::maxreal(cerr, fabs(c(3)+double(200)/double(120)));
    cerr = ap::maxreal(cerr, fabs(c(4)-double(25)/double(120)));
    cerr = ap::maxreal(cerr, fabs(c(5)+double(1)/double(120)));
    laguerrecoefficients(6, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-double(720)/double(720)));
    cerr = ap::maxreal(cerr, fabs(c(1)+double(4320)/double(720)));
    cerr = ap::maxreal(cerr, fabs(c(2)-double(5400)/double(720)));
    cerr = ap::maxreal(cerr, fabs(c(3)+double(2400)/double(720)));
    cerr = ap::maxreal(cerr, fabs(c(4)-double(450)/double(720)));
    cerr = ap::maxreal(cerr, fabs(c(5)+double(36)/double(720)));
    cerr = ap::maxreal(cerr, fabs(c(6)-double(1)/double(720)));
    
    //
    // Reporting
    //
    waserrors = ap::fp_greater(err,threshold)||ap::fp_greater(sumerr,threshold)||ap::fp_greater(cerr,threshold);
    if( !silent )
    {
        printf("TESTING CALCULATION OF THE LAGUERRE POLYNOMIALS\n");
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
bool testlaguerreunit_test_silent()
{
    bool result;

    result = testlaguerrecalculate(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testlaguerreunit_test()
{
    bool result;

    result = testlaguerrecalculate(false);
    return result;
}





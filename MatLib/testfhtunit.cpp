
 
#include <stdio.h>
#include "testfhtunit.h"

static void reffhtr1d(ap::real_1d_array& a, int n);
static void reffhtr1dinv(ap::real_1d_array& a, int n);

/*************************************************************************
Test
*************************************************************************/
bool testfht(bool silent)
{
    bool result;
    int n;
    int i;
    ap::real_1d_array r1;
    ap::real_1d_array r2;
    ap::real_1d_array r3;
    int maxn;
    double bidierr;
    double referr;
    double errtol;
    bool referrors;
    bool bidierrors;
    bool waserrors;

    maxn = 128;
    errtol = 100000*pow(double(maxn), double(3)/double(2))*ap::machineepsilon;
    bidierrors = false;
    referrors = false;
    waserrors = false;
    
    //
    // Test bi-directional error: norm(x-invFHT(FHT(x)))
    //
    bidierr = 0;
    for(n = 1; n <= maxn; n++)
    {
        
        //
        // FHT/invFHT
        //
        r1.setlength(n);
        r2.setlength(n);
        r3.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            r1(i) = 2*ap::randomreal()-1;
            r2(i) = r1(i);
            r3(i) = r1(i);
        }
        fhtr1d(r2, n);
        fhtr1dinv(r2, n);
        fhtr1dinv(r3, n);
        fhtr1d(r3, n);
        for(i = 0; i <= n-1; i++)
        {
            bidierr = ap::maxreal(bidierr, fabs(r1(i)-r2(i)));
            bidierr = ap::maxreal(bidierr, fabs(r1(i)-r3(i)));
        }
    }
    bidierrors = bidierrors||ap::fp_greater(bidierr,errtol);
    
    //
    // Test against reference O(N^2) implementation
    //
    referr = 0;
    for(n = 1; n <= maxn; n++)
    {
        
        //
        // FHT
        //
        r1.setlength(n);
        r2.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            r1(i) = 2*ap::randomreal()-1;
            r2(i) = r1(i);
        }
        fhtr1d(r1, n);
        reffhtr1d(r2, n);
        for(i = 0; i <= n-1; i++)
        {
            referr = ap::maxreal(referr, fabs(r1(i)-r2(i)));
        }
        
        //
        // inverse FHT
        //
        r1.setlength(n);
        r2.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            r1(i) = 2*ap::randomreal()-1;
            r2(i) = r1(i);
        }
        fhtr1dinv(r1, n);
        reffhtr1dinv(r2, n);
        for(i = 0; i <= n-1; i++)
        {
            referr = ap::maxreal(referr, fabs(r1(i)-r2(i)));
        }
    }
    referrors = referrors||ap::fp_greater(referr,errtol);
    
    //
    // end
    //
    waserrors = bidierrors||referrors;
    if( !silent )
    {
        printf("TESTING FHT\n");
        printf("FINAL RESULT:                             ");
        if( waserrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* BI-DIRECTIONAL TEST:                    ");
        if( bidierrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* AGAINST REFERENCE FHT:                  ");
        if( referrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        if( waserrors )
        {
            printf("TEST FAILED\n");
        }
        else
        {
            printf("TEST PASSED\n");
        }
    }
    result = !waserrors;
    return result;
}


/*************************************************************************
Reference FHT
*************************************************************************/
static void reffhtr1d(ap::real_1d_array& a, int n)
{
    ap::real_1d_array buf;
    int i;
    int j;
    double v;

    ap::ap_error::make_assertion(n>0, "RefFHTR1D: incorrect N!");
    buf.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        v = 0;
        for(j = 0; j <= n-1; j++)
        {
            v = v+a(j)*(cos(2*ap::pi()*i*j/n)+sin(2*ap::pi()*i*j/n));
        }
        buf(i) = v;
    }
    for(i = 0; i <= n-1; i++)
    {
        a(i) = buf(i);
    }
}


/*************************************************************************
Reference inverse FHT
*************************************************************************/
static void reffhtr1dinv(ap::real_1d_array& a, int n)
{
    int i;

    ap::ap_error::make_assertion(n>0, "RefFHTR1DInv: incorrect N!");
    reffhtr1d(a, n);
    for(i = 0; i <= n-1; i++)
    {
        a(i) = a(i)/n;
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testfhtunit_test_silent()
{
    bool result;

    result = testfht(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testfhtunit_test()
{
    bool result;

    result = testfht(false);
    return result;
}





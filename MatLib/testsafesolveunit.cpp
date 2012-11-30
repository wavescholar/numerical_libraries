
 
#include <stdio.h>
#include "testsafesolveunit.h"

static void rmatrixmakeacopy(const ap::real_2d_array& a,
     int m,
     int n,
     ap::real_2d_array& b);
static void cmatrixmakeacopy(const ap::complex_2d_array& a,
     int m,
     int n,
     ap::complex_2d_array& b);

/*************************************************************************
Main unittest subroutine
*************************************************************************/
bool testsafesolve(bool silent)
{
    bool result;
    int maxmn;
    double threshold;
    bool rerrors;
    bool cerrors;
    bool waserrors;
    bool isupper;
    int trans;
    bool isunit;
    double scalea;
    double growth;
    int i;
    int j;
    int n;
    int j1;
    int j2;
    ap::complex cv;
    ap::complex_2d_array ca;
    ap::complex_2d_array cea;
    ap::complex_2d_array ctmpa;
    ap::complex_1d_array cxs;
    ap::complex_1d_array cxe;
    double rv;
    ap::real_2d_array ra;
    ap::real_2d_array rea;
    ap::real_2d_array rtmpa;
    ap::real_1d_array rxs;
    ap::real_1d_array rxe;

    maxmn = 30;
    threshold = 100000*ap::machineepsilon;
    rerrors = false;
    cerrors = false;
    waserrors = false;
    
    //
    // Different problems: general tests
    //
    for(n = 1; n <= maxmn; n++)
    {
        
        //
        // test complex solver with well-conditioned matrix:
        // 1. generate A: fill off-diagonal elements with small values,
        //    diagonal elements are filled with larger values
        // 2. generate 'effective' A
        // 3. prepare task (exact X is stored in CXE, right part - in CXS),
        //    solve and compare CXS and CXE
        //
        isupper = ap::fp_greater(ap::randomreal(),0.5);
        trans = ap::randominteger(3);
        isunit = ap::fp_greater(ap::randomreal(),0.5);
        scalea = ap::randomreal()+0.5;
        ca.setlength(n, n);
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                if( i==j )
                {
                    ca(i,j).x = (2*ap::randominteger(2)-1)*(5+ap::randomreal());
                    ca(i,j).y = (2*ap::randominteger(2)-1)*(5+ap::randomreal());
                }
                else
                {
                    ca(i,j).x = 0.2*ap::randomreal()-0.1;
                    ca(i,j).y = 0.2*ap::randomreal()-0.1;
                }
            }
        }
        cmatrixmakeacopy(ca, n, n, ctmpa);
        for(i = 0; i <= n-1; i++)
        {
            if( isupper )
            {
                j1 = 0;
                j2 = i-1;
            }
            else
            {
                j1 = i+1;
                j2 = n-1;
            }
            for(j = j1; j <= j2; j++)
            {
                ctmpa(i,j) = 0;
            }
            if( isunit )
            {
                ctmpa(i,i) = 1;
            }
        }
        cea.setlength(n, n);
        for(i = 0; i <= n-1; i++)
        {
            if( trans==0 )
            {
                ap::vmove(&cea(i, 0), 1, &ctmpa(i, 0), 1, "N", ap::vlen(0,n-1), scalea);
            }
            if( trans==1 )
            {
                ap::vmove(&cea(0, i), cea.getstride(), &ctmpa(i, 0), 1, "N", ap::vlen(0,n-1), scalea);
            }
            if( trans==2 )
            {
                ap::vmove(&cea(0, i), cea.getstride(), &ctmpa(i, 0), 1, "Conj", ap::vlen(0,n-1), scalea);
            }
        }
        cxe.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            cxe(i).x = 2*ap::randomreal()-1;
            cxe(i).y = 2*ap::randomreal()-1;
        }
        cxs.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            cv = ap::vdotproduct(&cea(i, 0), 1, "N", &cxe(0), 1, "N", ap::vlen(0,n-1));
            cxs(i) = cv;
        }
        if( cmatrixscaledtrsafesolve(ca, scalea, n, cxs, isupper, trans, isunit, sqrt(ap::maxrealnumber)) )
        {
            for(i = 0; i <= n-1; i++)
            {
                cerrors = cerrors||ap::fp_greater(ap::abscomplex(cxs(i)-cxe(i)),threshold);
            }
        }
        else
        {
            cerrors = true;
        }
        
        //
        // same with real
        //
        isupper = ap::fp_greater(ap::randomreal(),0.5);
        trans = ap::randominteger(2);
        isunit = ap::fp_greater(ap::randomreal(),0.5);
        scalea = ap::randomreal()+0.5;
        ra.setlength(n, n);
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                if( i==j )
                {
                    ra(i,j) = (2*ap::randominteger(2)-1)*(5+ap::randomreal());
                }
                else
                {
                    ra(i,j) = 0.2*ap::randomreal()-0.1;
                }
            }
        }
        rmatrixmakeacopy(ra, n, n, rtmpa);
        for(i = 0; i <= n-1; i++)
        {
            if( isupper )
            {
                j1 = 0;
                j2 = i-1;
            }
            else
            {
                j1 = i+1;
                j2 = n-1;
            }
            for(j = j1; j <= j2; j++)
            {
                rtmpa(i,j) = 0;
            }
            if( isunit )
            {
                rtmpa(i,i) = 1;
            }
        }
        rea.setlength(n, n);
        for(i = 0; i <= n-1; i++)
        {
            if( trans==0 )
            {
                ap::vmove(&rea(i, 0), 1, &rtmpa(i, 0), 1, ap::vlen(0,n-1), scalea);
            }
            if( trans==1 )
            {
                ap::vmove(&rea(0, i), rea.getstride(), &rtmpa(i, 0), 1, ap::vlen(0,n-1), scalea);
            }
        }
        rxe.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            rxe(i) = 2*ap::randomreal()-1;
        }
        rxs.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            rv = ap::vdotproduct(&rea(i, 0), 1, &rxe(0), 1, ap::vlen(0,n-1));
            rxs(i) = rv;
        }
        if( rmatrixscaledtrsafesolve(ra, scalea, n, rxs, isupper, trans, isunit, sqrt(ap::maxrealnumber)) )
        {
            for(i = 0; i <= n-1; i++)
            {
                rerrors = rerrors||ap::fp_greater(fabs(rxs(i)-rxe(i)),threshold);
            }
        }
        else
        {
            rerrors = true;
        }
    }
    
    //
    // Special test with diagonal ill-conditioned matrix:
    // * ability to solve it when resulting growth is less than threshold
    // * ability to stop solve when resulting growth is greater than threshold
    //
    // A = diag(1, 1/growth)
    // b = (1, 0.5)
    //
    n = 2;
    growth = 10;
    ca.setlength(n, n);
    ca(0,0) = 1;
    ca(0,1) = 0;
    ca(1,0) = 0;
    ca(1,1) = 1/growth;
    cxs.setlength(n);
    cxs(0) = 1.0;
    cxs(1) = 0.5;
    cerrors = cerrors||!cmatrixscaledtrsafesolve(ca, 1.0, n, cxs, ap::fp_greater(ap::randomreal(),0.5), ap::randominteger(3), false, 1.05*ap::maxreal(ap::abscomplex(cxs(1))*growth, 1.0));
    cerrors = cerrors||!cmatrixscaledtrsafesolve(ca, 1.0, n, cxs, ap::fp_greater(ap::randomreal(),0.5), ap::randominteger(3), false, 0.95*ap::maxreal(ap::abscomplex(cxs(1))*growth, 1.0));
    ra.setlength(n, n);
    ra(0,0) = 1;
    ra(0,1) = 0;
    ra(1,0) = 0;
    ra(1,1) = 1/growth;
    rxs.setlength(n);
    rxs(0) = 1.0;
    rxs(1) = 0.5;
    rerrors = rerrors||!rmatrixscaledtrsafesolve(ra, 1.0, n, rxs, ap::fp_greater(ap::randomreal(),0.5), ap::randominteger(2), false, 1.05*ap::maxreal(fabs(rxs(1))*growth, 1.0));
    rerrors = rerrors||!rmatrixscaledtrsafesolve(ra, 1.0, n, rxs, ap::fp_greater(ap::randomreal(),0.5), ap::randominteger(2), false, 0.95*ap::maxreal(fabs(rxs(1))*growth, 1.0));
    
    //
    // Special test with diagonal degenerate matrix:
    // * ability to solve it when resulting growth is less than threshold
    // * ability to stop solve when resulting growth is greater than threshold
    //
    // A = diag(1, 0)
    // b = (1, 0.5)
    //
    n = 2;
    ca.setlength(n, n);
    ca(0,0) = 1;
    ca(0,1) = 0;
    ca(1,0) = 0;
    ca(1,1) = 0;
    cxs.setlength(n);
    cxs(0) = 1.0;
    cxs(1) = 0.5;
    cerrors = cerrors||cmatrixscaledtrsafesolve(ca, 1.0, n, cxs, ap::fp_greater(ap::randomreal(),0.5), ap::randominteger(3), false, sqrt(ap::maxrealnumber));
    ra.setlength(n, n);
    ra(0,0) = 1;
    ra(0,1) = 0;
    ra(1,0) = 0;
    ra(1,1) = 0;
    rxs.setlength(n);
    rxs(0) = 1.0;
    rxs(1) = 0.5;
    rerrors = rerrors||rmatrixscaledtrsafesolve(ra, 1.0, n, rxs, ap::fp_greater(ap::randomreal(),0.5), ap::randominteger(2), false, sqrt(ap::maxrealnumber));
    
    //
    // report
    //
    waserrors = rerrors||cerrors;
    if( !silent )
    {
        printf("TESTING SAFE TR SOLVER\n");
        printf("REAL:                                    ");
        if( !rerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("COMPLEX:                                 ");
        if( !cerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        if( waserrors )
        {
            printf("TEST FAILED\n");
        }
        else
        {
            printf("TEST PASSED\n");
        }
        printf("\n\n");
    }
    result = !waserrors;
    return result;
}


/*************************************************************************
Copy
*************************************************************************/
static void rmatrixmakeacopy(const ap::real_2d_array& a,
     int m,
     int n,
     ap::real_2d_array& b)
{
    int i;
    int j;

    b.setbounds(0, m-1, 0, n-1);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            b(i,j) = a(i,j);
        }
    }
}


/*************************************************************************
Copy
*************************************************************************/
static void cmatrixmakeacopy(const ap::complex_2d_array& a,
     int m,
     int n,
     ap::complex_2d_array& b)
{
    int i;
    int j;

    b.setbounds(0, m-1, 0, n-1);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            b(i,j) = a(i,j);
        }
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testsafesolveunit_test_silent()
{
    bool result;

    result = testsafesolve(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testsafesolveunit_test()
{
    bool result;

    result = testsafesolve(false);
    return result;
}





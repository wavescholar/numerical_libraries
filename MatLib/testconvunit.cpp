
 
#include <stdio.h>
#include "testconvunit.h"

static void refconvc1d(const ap::complex_1d_array& a,
     int m,
     const ap::complex_1d_array& b,
     int n,
     ap::complex_1d_array& r);
static void refconvc1dcircular(const ap::complex_1d_array& a,
     int m,
     const ap::complex_1d_array& b,
     int n,
     ap::complex_1d_array& r);
static void refconvr1d(const ap::real_1d_array& a,
     int m,
     const ap::real_1d_array& b,
     int n,
     ap::real_1d_array& r);
static void refconvr1dcircular(const ap::real_1d_array& a,
     int m,
     const ap::real_1d_array& b,
     int n,
     ap::real_1d_array& r);

/*************************************************************************
Test
*************************************************************************/
bool testconv(bool silent)
{
    bool result;
    int m;
    int n;
    int i;
    int rkind;
    int circkind;
    ap::real_1d_array ra;
    ap::real_1d_array rb;
    ap::real_1d_array rr1;
    ap::real_1d_array rr2;
    ap::complex_1d_array ca;
    ap::complex_1d_array cb;
    ap::complex_1d_array cr1;
    ap::complex_1d_array cr2;
    int maxn;
    double referr;
    double refrerr;
    double inverr;
    double invrerr;
    double errtol;
    bool referrors;
    bool refrerrors;
    bool inverrors;
    bool invrerrors;
    bool waserrors;

    maxn = 32;
    errtol = 100000*pow(double(maxn), double(3)/double(2))*ap::machineepsilon;
    referrors = false;
    refrerrors = false;
    inverrors = false;
    invrerrors = false;
    waserrors = false;
    
    //
    // Test against reference O(N^2) implementation.
    //
    // Automatic ConvC1D() and different algorithms of ConvC1DX() are tested.
    //
    referr = 0;
    refrerr = 0;
    for(m = 1; m <= maxn; m++)
    {
        for(n = 1; n <= maxn; n++)
        {
            for(circkind = 0; circkind <= 1; circkind++)
            {
                for(rkind = -3; rkind <= 1; rkind++)
                {
                    
                    //
                    // skip impossible combinations of parameters:
                    // * circular convolution, M<N, RKind<>-3 - internal subroutine does not support M<N.
                    //
                    if( circkind!=0&&m<n&&rkind!=-3 )
                    {
                        continue;
                    }
                    
                    //
                    // Complex convolution
                    //
                    ca.setlength(m);
                    for(i = 0; i <= m-1; i++)
                    {
                        ca(i).x = 2*ap::randomreal()-1;
                        ca(i).y = 2*ap::randomreal()-1;
                    }
                    cb.setlength(n);
                    for(i = 0; i <= n-1; i++)
                    {
                        cb(i).x = 2*ap::randomreal()-1;
                        cb(i).y = 2*ap::randomreal()-1;
                    }
                    cr1.setlength(1);
                    if( rkind==-3 )
                    {
                        
                        //
                        // test wrapper subroutine:
                        // * circular/non-circular
                        //
                        if( circkind==0 )
                        {
                            convc1d(ca, m, cb, n, cr1);
                        }
                        else
                        {
                            convc1dcircular(ca, m, cb, n, cr1);
                        }
                    }
                    else
                    {
                        
                        //
                        // test internal subroutine
                        //
                        if( m>=n )
                        {
                            
                            //
                            // test internal subroutine:
                            // * circular/non-circular mode
                            //
                            convc1dx(ca, m, cb, n, circkind!=0, rkind, 0, cr1);
                        }
                        else
                        {
                            
                            //
                            // test internal subroutine - circular mode only
                            //
                            ap::ap_error::make_assertion(circkind==0, "Convolution test: internal error!");
                            convc1dx(cb, n, ca, m, false, rkind, 0, cr1);
                        }
                    }
                    if( circkind==0 )
                    {
                        refconvc1d(ca, m, cb, n, cr2);
                    }
                    else
                    {
                        refconvc1dcircular(ca, m, cb, n, cr2);
                    }
                    if( circkind==0 )
                    {
                        for(i = 0; i <= m+n-2; i++)
                        {
                            referr = ap::maxreal(referr, ap::abscomplex(cr1(i)-cr2(i)));
                        }
                    }
                    else
                    {
                        for(i = 0; i <= m-1; i++)
                        {
                            referr = ap::maxreal(referr, ap::abscomplex(cr1(i)-cr2(i)));
                        }
                    }
                    
                    //
                    // Real convolution
                    //
                    ra.setlength(m);
                    for(i = 0; i <= m-1; i++)
                    {
                        ra(i) = 2*ap::randomreal()-1;
                    }
                    rb.setlength(n);
                    for(i = 0; i <= n-1; i++)
                    {
                        rb(i) = 2*ap::randomreal()-1;
                    }
                    rr1.setlength(1);
                    if( rkind==-3 )
                    {
                        
                        //
                        // test wrapper subroutine:
                        // * circular/non-circular
                        //
                        if( circkind==0 )
                        {
                            convr1d(ra, m, rb, n, rr1);
                        }
                        else
                        {
                            convr1dcircular(ra, m, rb, n, rr1);
                        }
                    }
                    else
                    {
                        if( m>=n )
                        {
                            
                            //
                            // test internal subroutine:
                            // * circular/non-circular mode
                            //
                            convr1dx(ra, m, rb, n, circkind!=0, rkind, 0, rr1);
                        }
                        else
                        {
                            
                            //
                            // test internal subroutine - non-circular mode only
                            //
                            convr1dx(rb, n, ra, m, circkind!=0, rkind, 0, rr1);
                        }
                    }
                    if( circkind==0 )
                    {
                        refconvr1d(ra, m, rb, n, rr2);
                    }
                    else
                    {
                        refconvr1dcircular(ra, m, rb, n, rr2);
                    }
                    if( circkind==0 )
                    {
                        for(i = 0; i <= m+n-2; i++)
                        {
                            refrerr = ap::maxreal(refrerr, fabs(rr1(i)-rr2(i)));
                        }
                    }
                    else
                    {
                        for(i = 0; i <= m-1; i++)
                        {
                            refrerr = ap::maxreal(refrerr, fabs(rr1(i)-rr2(i)));
                        }
                    }
                }
            }
        }
    }
    referrors = referrors||ap::fp_greater(referr,errtol);
    refrerrors = refrerrors||ap::fp_greater(refrerr,errtol);
    
    //
    // Test inverse convolution
    //
    inverr = 0;
    invrerr = 0;
    for(m = 1; m <= maxn; m++)
    {
        for(n = 1; n <= maxn; n++)
        {
            
            //
            // Complex circilar and non-circular
            //
            ca.setlength(m);
            for(i = 0; i <= m-1; i++)
            {
                ca(i).x = 2*ap::randomreal()-1;
                ca(i).y = 2*ap::randomreal()-1;
            }
            cb.setlength(n);
            for(i = 0; i <= n-1; i++)
            {
                cb(i).x = 2*ap::randomreal()-1;
                cb(i).y = 2*ap::randomreal()-1;
            }
            cr1.setlength(1);
            cr2.setlength(1);
            convc1d(ca, m, cb, n, cr2);
            convc1dinv(cr2, m+n-1, cb, n, cr1);
            for(i = 0; i <= m-1; i++)
            {
                inverr = ap::maxreal(inverr, ap::abscomplex(cr1(i)-ca(i)));
            }
            cr1.setlength(1);
            cr2.setlength(1);
            convc1dcircular(ca, m, cb, n, cr2);
            convc1dcircularinv(cr2, m, cb, n, cr1);
            for(i = 0; i <= m-1; i++)
            {
                inverr = ap::maxreal(inverr, ap::abscomplex(cr1(i)-ca(i)));
            }
            
            //
            // Real circilar and non-circular
            //
            ra.setlength(m);
            for(i = 0; i <= m-1; i++)
            {
                ra(i) = 2*ap::randomreal()-1;
            }
            rb.setlength(n);
            for(i = 0; i <= n-1; i++)
            {
                rb(i) = 2*ap::randomreal()-1;
            }
            rr1.setlength(1);
            rr2.setlength(1);
            convr1d(ra, m, rb, n, rr2);
            convr1dinv(rr2, m+n-1, rb, n, rr1);
            for(i = 0; i <= m-1; i++)
            {
                invrerr = ap::maxreal(invrerr, fabs(rr1(i)-ra(i)));
            }
            rr1.setlength(1);
            rr2.setlength(1);
            convr1dcircular(ra, m, rb, n, rr2);
            convr1dcircularinv(rr2, m, rb, n, rr1);
            for(i = 0; i <= m-1; i++)
            {
                invrerr = ap::maxreal(invrerr, fabs(rr1(i)-ra(i)));
            }
        }
    }
    inverrors = inverrors||ap::fp_greater(inverr,errtol);
    invrerrors = invrerrors||ap::fp_greater(invrerr,errtol);
    
    //
    // end
    //
    waserrors = referrors||refrerrors||inverrors||invrerrors;
    if( !silent )
    {
        printf("TESTING CONVOLUTION\n");
        printf("FINAL RESULT:                             ");
        if( waserrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* AGAINST REFERENCE COMPLEX CONV:         ");
        if( referrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* AGAINST REFERENCE REAL CONV:            ");
        if( refrerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* COMPLEX INVERSE:                        ");
        if( inverrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* REAL INVERSE:                           ");
        if( invrerrors )
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
Reference implementation
*************************************************************************/
static void refconvc1d(const ap::complex_1d_array& a,
     int m,
     const ap::complex_1d_array& b,
     int n,
     ap::complex_1d_array& r)
{
    int i;
    ap::complex v;

    r.setlength(m+n-1);
    for(i = 0; i <= m+n-2; i++)
    {
        r(i) = 0;
    }
    for(i = 0; i <= m-1; i++)
    {
        v = a(i);
        ap::vadd(&r(i), 1, &b(0), 1, "N", ap::vlen(i,i+n-1), v);
    }
}


/*************************************************************************
Reference implementation
*************************************************************************/
static void refconvc1dcircular(const ap::complex_1d_array& a,
     int m,
     const ap::complex_1d_array& b,
     int n,
     ap::complex_1d_array& r)
{
    int i1;
    int i2;
    int j2;
    ap::complex_1d_array buf;

    refconvc1d(a, m, b, n, buf);
    r.setlength(m);
    ap::vmove(&r(0), 1, &buf(0), 1, "N", ap::vlen(0,m-1));
    i1 = m;
    while(i1<=m+n-2)
    {
        i2 = ap::minint(i1+m-1, m+n-2);
        j2 = i2-i1;
        ap::vadd(&r(0), 1, &buf(i1), 1, "N", ap::vlen(0,j2));
        i1 = i1+m;
    }
}


/*************************************************************************
Reference FFT
*************************************************************************/
static void refconvr1d(const ap::real_1d_array& a,
     int m,
     const ap::real_1d_array& b,
     int n,
     ap::real_1d_array& r)
{
    int i;
    double v;

    r.setlength(m+n-1);
    for(i = 0; i <= m+n-2; i++)
    {
        r(i) = 0;
    }
    for(i = 0; i <= m-1; i++)
    {
        v = a(i);
        ap::vadd(&r(i), 1, &b(0), 1, ap::vlen(i,i+n-1), v);
    }
}


/*************************************************************************
Reference implementation
*************************************************************************/
static void refconvr1dcircular(const ap::real_1d_array& a,
     int m,
     const ap::real_1d_array& b,
     int n,
     ap::real_1d_array& r)
{
    int i1;
    int i2;
    int j2;
    ap::real_1d_array buf;

    refconvr1d(a, m, b, n, buf);
    r.setlength(m);
    ap::vmove(&r(0), 1, &buf(0), 1, ap::vlen(0,m-1));
    i1 = m;
    while(i1<=m+n-2)
    {
        i2 = ap::minint(i1+m-1, m+n-2);
        j2 = i2-i1;
        ap::vadd(&r(0), 1, &buf(i1), 1, ap::vlen(0,j2));
        i1 = i1+m;
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testconvunit_test_silent()
{
    bool result;

    result = testconv(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testconvunit_test()
{
    bool result;

    result = testconv(false);
    return result;
}





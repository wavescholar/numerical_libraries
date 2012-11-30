
 
#include <stdio.h>
#include "testcorrunit.h"

static void refcorrc1d(const ap::complex_1d_array& signal,
     int n,
     const ap::complex_1d_array& pattern,
     int m,
     ap::complex_1d_array& r);
static void refcorrc1dcircular(const ap::complex_1d_array& signal,
     int n,
     const ap::complex_1d_array& pattern,
     int m,
     ap::complex_1d_array& r);
static void refcorrr1d(const ap::real_1d_array& signal,
     int n,
     const ap::real_1d_array& pattern,
     int m,
     ap::real_1d_array& r);
static void refcorrr1dcircular(const ap::real_1d_array& signal,
     int n,
     const ap::real_1d_array& pattern,
     int m,
     ap::real_1d_array& r);
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
bool testcorr(bool silent)
{
    bool result;
    int m;
    int n;
    int i;
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
    referr = 0;
    refrerr = 0;
    for(m = 1; m <= maxn; m++)
    {
        for(n = 1; n <= maxn; n++)
        {
            
            //
            // Complex correlation
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
            corrc1d(ca, m, cb, n, cr1);
            refcorrc1d(ca, m, cb, n, cr2);
            for(i = 0; i <= m+n-2; i++)
            {
                referr = ap::maxreal(referr, ap::abscomplex(cr1(i)-cr2(i)));
            }
            cr1.setlength(1);
            corrc1dcircular(ca, m, cb, n, cr1);
            refcorrc1dcircular(ca, m, cb, n, cr2);
            for(i = 0; i <= m-1; i++)
            {
                referr = ap::maxreal(referr, ap::abscomplex(cr1(i)-cr2(i)));
            }
            
            //
            // Real correlation
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
            corrr1d(ra, m, rb, n, rr1);
            refcorrr1d(ra, m, rb, n, rr2);
            for(i = 0; i <= m+n-2; i++)
            {
                refrerr = ap::maxreal(refrerr, fabs(rr1(i)-rr2(i)));
            }
            rr1.setlength(1);
            corrr1dcircular(ra, m, rb, n, rr1);
            refcorrr1dcircular(ra, m, rb, n, rr2);
            for(i = 0; i <= m-1; i++)
            {
                refrerr = ap::maxreal(refrerr, fabs(rr1(i)-rr2(i)));
            }
        }
    }
    referrors = referrors||ap::fp_greater(referr,errtol);
    refrerrors = refrerrors||ap::fp_greater(refrerr,errtol);
    
    //
    // end
    //
    waserrors = referrors||refrerrors;
    if( !silent )
    {
        printf("TESTING CORRELATION\n");
        printf("FINAL RESULT:                             ");
        if( waserrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* AGAINST REFERENCE COMPLEX CORR:         ");
        if( referrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* AGAINST REFERENCE REAL CORR:            ");
        if( refrerrors )
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
static void refcorrc1d(const ap::complex_1d_array& signal,
     int n,
     const ap::complex_1d_array& pattern,
     int m,
     ap::complex_1d_array& r)
{
    int i;
    int j;
    ap::complex v;
    ap::complex_1d_array s;

    s.setlength(m+n-1);
    ap::vmove(&s(0), 1, &signal(0), 1, "N", ap::vlen(0,n-1));
    for(i = n; i <= m+n-2; i++)
    {
        s(i) = 0;
    }
    r.setlength(m+n-1);
    for(i = 0; i <= n-1; i++)
    {
        v = 0;
        for(j = 0; j <= m-1; j++)
        {
            if( i+j>=n )
            {
                break;
            }
            v = v+ap::conj(pattern(j))*s(i+j);
        }
        r(i) = v;
    }
    for(i = 1; i <= m-1; i++)
    {
        v = 0;
        for(j = i; j <= m-1; j++)
        {
            v = v+ap::conj(pattern(j))*s(j-i);
        }
        r(m+n-1-i) = v;
    }
}


/*************************************************************************
Reference implementation
*************************************************************************/
static void refcorrc1dcircular(const ap::complex_1d_array& signal,
     int n,
     const ap::complex_1d_array& pattern,
     int m,
     ap::complex_1d_array& r)
{
    int i;
    int j;
    ap::complex v;

    r.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        v = 0;
        for(j = 0; j <= m-1; j++)
        {
            v = v+ap::conj(pattern(j))*signal((i+j)%n);
        }
        r(i) = v;
    }
}


/*************************************************************************
Reference implementation
*************************************************************************/
static void refcorrr1d(const ap::real_1d_array& signal,
     int n,
     const ap::real_1d_array& pattern,
     int m,
     ap::real_1d_array& r)
{
    int i;
    int j;
    double v;
    ap::real_1d_array s;

    s.setlength(m+n-1);
    ap::vmove(&s(0), 1, &signal(0), 1, ap::vlen(0,n-1));
    for(i = n; i <= m+n-2; i++)
    {
        s(i) = 0;
    }
    r.setlength(m+n-1);
    for(i = 0; i <= n-1; i++)
    {
        v = 0;
        for(j = 0; j <= m-1; j++)
        {
            if( i+j>=n )
            {
                break;
            }
            v = v+pattern(j)*s(i+j);
        }
        r(i) = v;
    }
    for(i = 1; i <= m-1; i++)
    {
        v = 0;
        for(j = i; j <= m-1; j++)
        {
            v = v+pattern(j)*s(-i+j);
        }
        r(m+n-1-i) = v;
    }
}


/*************************************************************************
Reference implementation
*************************************************************************/
static void refcorrr1dcircular(const ap::real_1d_array& signal,
     int n,
     const ap::real_1d_array& pattern,
     int m,
     ap::real_1d_array& r)
{
    int i;
    int j;
    double v;

    r.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        v = 0;
        for(j = 0; j <= m-1; j++)
        {
            v = v+pattern(j)*signal((i+j)%n);
        }
        r(i) = v;
    }
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
bool testcorrunit_test_silent()
{
    bool result;

    result = testcorr(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testcorrunit_test()
{
    bool result;

    result = testcorr(false);
    return result;
}





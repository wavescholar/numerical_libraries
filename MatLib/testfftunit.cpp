
 
#include <stdio.h>
#include "testfftunit.h"

static void reffftc1d(ap::complex_1d_array& a, int n);
static void reffftc1dinv(ap::complex_1d_array& a, int n);
static void refinternalcfft(ap::real_1d_array& a, int nn, bool inversefft);
static void refinternalrfft(const ap::real_1d_array& a,
     int nn,
     ap::complex_1d_array& f);

/*************************************************************************
Test
*************************************************************************/
bool testfft(bool silent)
{
    bool result;
    int n;
    int i;
    int k;
    ap::complex_1d_array a1;
    ap::complex_1d_array a2;
    ap::complex_1d_array a3;
    ap::real_1d_array r1;
    ap::real_1d_array r2;
    ap::real_1d_array buf;
    ftplan plan;
    int maxn;
    double bidierr;
    double bidirerr;
    double referr;
    double refrerr;
    double reinterr;
    double errtol;
    bool referrors;
    bool bidierrors;
    bool refrerrors;
    bool bidirerrors;
    bool reinterrors;
    bool waserrors;

    maxn = 128;
    errtol = 100000*pow(double(maxn), double(3)/double(2))*ap::machineepsilon;
    bidierrors = false;
    referrors = false;
    bidirerrors = false;
    refrerrors = false;
    reinterrors = false;
    waserrors = false;
    
    //
    // Test bi-directional error: norm(x-invFFT(FFT(x)))
    //
    bidierr = 0;
    bidirerr = 0;
    for(n = 1; n <= maxn; n++)
    {
        
        //
        // Complex FFT/invFFT
        //
        a1.setlength(n);
        a2.setlength(n);
        a3.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            a1(i).x = 2*ap::randomreal()-1;
            a1(i).y = 2*ap::randomreal()-1;
            a2(i) = a1(i);
            a3(i) = a1(i);
        }
        fftc1d(a2, n);
        fftc1dinv(a2, n);
        fftc1dinv(a3, n);
        fftc1d(a3, n);
        for(i = 0; i <= n-1; i++)
        {
            bidierr = ap::maxreal(bidierr, ap::abscomplex(a1(i)-a2(i)));
            bidierr = ap::maxreal(bidierr, ap::abscomplex(a1(i)-a3(i)));
        }
        
        //
        // Real
        //
        r1.setlength(n);
        r2.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            r1(i) = 2*ap::randomreal()-1;
            r2(i) = r1(i);
        }
        fftr1d(r2, n, a1);
        ap::vmul(&r2(0), 1, ap::vlen(0,n-1), 0);
        fftr1dinv(a1, n, r2);
        for(i = 0; i <= n-1; i++)
        {
            bidirerr = ap::maxreal(bidirerr, ap::abscomplex(r1(i)-r2(i)));
        }
    }
    bidierrors = bidierrors||ap::fp_greater(bidierr,errtol);
    bidirerrors = bidirerrors||ap::fp_greater(bidirerr,errtol);
    
    //
    // Test against reference O(N^2) implementation
    //
    referr = 0;
    refrerr = 0;
    for(n = 1; n <= maxn; n++)
    {
        
        //
        // Complex FFT
        //
        a1.setlength(n);
        a2.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            a1(i).x = 2*ap::randomreal()-1;
            a1(i).y = 2*ap::randomreal()-1;
            a2(i) = a1(i);
        }
        fftc1d(a1, n);
        reffftc1d(a2, n);
        for(i = 0; i <= n-1; i++)
        {
            referr = ap::maxreal(referr, ap::abscomplex(a1(i)-a2(i)));
        }
        
        //
        // Complex inverse FFT
        //
        a1.setlength(n);
        a2.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            a1(i).x = 2*ap::randomreal()-1;
            a1(i).y = 2*ap::randomreal()-1;
            a2(i) = a1(i);
        }
        fftc1dinv(a1, n);
        reffftc1dinv(a2, n);
        for(i = 0; i <= n-1; i++)
        {
            referr = ap::maxreal(referr, ap::abscomplex(a1(i)-a2(i)));
        }
        
        //
        // Real forward/inverse FFT:
        // * calculate and check forward FFT
        // * use precalculated FFT to check backward FFT
        //   fill unused parts of frequencies array with random numbers
        //   to ensure that they are not really used
        //
        r1.setlength(n);
        r2.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            r1(i) = 2*ap::randomreal()-1;
            r2(i) = r1(i);
        }
        fftr1d(r1, n, a1);
        refinternalrfft(r2, n, a2);
        for(i = 0; i <= n-1; i++)
        {
            refrerr = ap::maxreal(refrerr, ap::abscomplex(a1(i)-a2(i)));
        }
        a3.setlength(ap::ifloor(double(n)/double(2))+1);
        for(i = 0; i <= ap::ifloor(double(n)/double(2)); i++)
        {
            a3(i) = a2(i);
        }
        a3(0).y = 2*ap::randomreal()-1;
        if( n%2==0 )
        {
            a3(ap::ifloor(double(n)/double(2))).y = 2*ap::randomreal()-1;
        }
        for(i = 0; i <= n-1; i++)
        {
            r1(i) = 0;
        }
        fftr1dinv(a3, n, r1);
        for(i = 0; i <= n-1; i++)
        {
            refrerr = ap::maxreal(refrerr, fabs(r2(i)-r1(i)));
        }
    }
    referrors = referrors||ap::fp_greater(referr,errtol);
    refrerrors = refrerrors||ap::fp_greater(refrerr,errtol);
    
    //
    // test internal real even FFT
    //
    reinterr = 0;
    for(k = 1; k <= maxn/2; k++)
    {
        n = 2*k;
        
        //
        // Real forward FFT
        //
        r1.setlength(n);
        r2.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            r1(i) = 2*ap::randomreal()-1;
            r2(i) = r1(i);
        }
        ftbasegeneratecomplexfftplan(n/2, plan);
        buf.setlength(n);
        fftr1dinternaleven(r1, n, buf, plan);
        refinternalrfft(r2, n, a2);
        reinterr = ap::maxreal(reinterr, fabs(r1(0)-a2(0).x));
        reinterr = ap::maxreal(reinterr, fabs(r1(1)-a2(n/2).x));
        for(i = 1; i <= n/2-1; i++)
        {
            reinterr = ap::maxreal(reinterr, fabs(r1(2*i+0)-a2(i).x));
            reinterr = ap::maxreal(reinterr, fabs(r1(2*i+1)-a2(i).y));
        }
        
        //
        // Real backward FFT
        //
        r1.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            r1(i) = 2*ap::randomreal()-1;
        }
        a2.setlength(ap::ifloor(double(n)/double(2))+1);
        a2(0) = r1(0);
        for(i = 1; i <= ap::ifloor(double(n)/double(2))-1; i++)
        {
            a2(i).x = r1(2*i+0);
            a2(i).y = r1(2*i+1);
        }
        a2(ap::ifloor(double(n)/double(2))) = r1(1);
        ftbasegeneratecomplexfftplan(n/2, plan);
        buf.setlength(n);
        fftr1dinvinternaleven(r1, n, buf, plan);
        fftr1dinv(a2, n, r2);
        for(i = 0; i <= n-1; i++)
        {
            reinterr = ap::maxreal(reinterr, fabs(r1(i)-r2(i)));
        }
    }
    reinterrors = reinterrors||ap::fp_greater(reinterr,errtol);
    
    //
    // end
    //
    waserrors = bidierrors||bidirerrors||referrors||refrerrors||reinterrors;
    if( !silent )
    {
        printf("TESTING FFT\n");
        printf("FINAL RESULT:                             ");
        if( waserrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* BI-DIRECTIONAL COMPLEX TEST:            ");
        if( bidierrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* AGAINST REFERENCE COMPLEX FFT:          ");
        if( referrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* BI-DIRECTIONAL REAL TEST:               ");
        if( bidirerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* AGAINST REFERENCE REAL FFT:             ");
        if( refrerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* INTERNAL EVEN FFT:                      ");
        if( reinterrors )
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
Reference FFT
*************************************************************************/
static void reffftc1d(ap::complex_1d_array& a, int n)
{
    ap::real_1d_array buf;
    int i;

    ap::ap_error::make_assertion(n>0, "FFTC1D: incorrect N!");
    buf.setlength(2*n);
    for(i = 0; i <= n-1; i++)
    {
        buf(2*i+0) = a(i).x;
        buf(2*i+1) = a(i).y;
    }
    refinternalcfft(buf, n, false);
    for(i = 0; i <= n-1; i++)
    {
        a(i).x = buf(2*i+0);
        a(i).y = buf(2*i+1);
    }
}


/*************************************************************************
Reference inverse FFT
*************************************************************************/
static void reffftc1dinv(ap::complex_1d_array& a, int n)
{
    ap::real_1d_array buf;
    int i;

    ap::ap_error::make_assertion(n>0, "FFTC1DInv: incorrect N!");
    buf.setlength(2*n);
    for(i = 0; i <= n-1; i++)
    {
        buf(2*i+0) = a(i).x;
        buf(2*i+1) = a(i).y;
    }
    refinternalcfft(buf, n, true);
    for(i = 0; i <= n-1; i++)
    {
        a(i).x = buf(2*i+0);
        a(i).y = buf(2*i+1);
    }
}


/*************************************************************************
Internal complex FFT stub.
Uses straightforward formula with O(N^2) complexity.
*************************************************************************/
static void refinternalcfft(ap::real_1d_array& a, int nn, bool inversefft)
{
    ap::real_1d_array tmp;
    int i;
    int j;
    int k;
    double hre;
    double him;
    double c;
    double s;
    double re;
    double im;

    tmp.setbounds(0, 2*nn-1);
    if( !inversefft )
    {
        for(i = 0; i <= nn-1; i++)
        {
            hre = 0;
            him = 0;
            for(k = 0; k <= nn-1; k++)
            {
                re = a(2*k);
                im = a(2*k+1);
                c = cos(-2*ap::pi()*k*i/nn);
                s = sin(-2*ap::pi()*k*i/nn);
                hre = hre+c*re-s*im;
                him = him+c*im+s*re;
            }
            tmp(2*i) = hre;
            tmp(2*i+1) = him;
        }
        for(i = 0; i <= 2*nn-1; i++)
        {
            a(i) = tmp(i);
        }
    }
    else
    {
        for(k = 0; k <= nn-1; k++)
        {
            hre = 0;
            him = 0;
            for(i = 0; i <= nn-1; i++)
            {
                re = a(2*i);
                im = a(2*i+1);
                c = cos(2*ap::pi()*k*i/nn);
                s = sin(2*ap::pi()*k*i/nn);
                hre = hre+c*re-s*im;
                him = him+c*im+s*re;
            }
            tmp(2*k) = hre/nn;
            tmp(2*k+1) = him/nn;
        }
        for(i = 0; i <= 2*nn-1; i++)
        {
            a(i) = tmp(i);
        }
    }
}


/*************************************************************************
Internal real FFT stub.
Uses straightforward formula with O(N^2) complexity.
*************************************************************************/
static void refinternalrfft(const ap::real_1d_array& a,
     int nn,
     ap::complex_1d_array& f)
{
    ap::real_1d_array tmp;
    int i;

    tmp.setbounds(0, 2*nn-1);
    for(i = 0; i <= nn-1; i++)
    {
        tmp(2*i) = a(i);
        tmp(2*i+1) = 0;
    }
    refinternalcfft(tmp, nn, false);
    f.setlength(nn);
    for(i = 0; i <= nn-1; i++)
    {
        f(i).x = tmp(2*i+0);
        f(i).y = tmp(2*i+1);
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testfftunit_test_silent()
{
    bool result;

    result = testfft(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testfftunit_test()
{
    bool result;

    result = testfft(false);
    return result;
}





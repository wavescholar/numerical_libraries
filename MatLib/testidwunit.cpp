
 
#include <stdio.h>
#include "testidwunit.h"

static void unset2d(ap::complex_2d_array& a);
static void unset1d(ap::real_1d_array& a);
static void testxy(const ap::real_2d_array& xy,
     const int& n,
     const int& nx,
     const int& d,
     const int& nq,
     const int& nw,
     bool& idwerrors);
static void testrxy(const ap::real_2d_array& xy,
     const int& n,
     const int& nx,
     const double& r,
     bool& idwerrors);
static void testdegree(const int& n,
     const int& nx,
     const int& d,
     const int& dtask,
     bool& idwerrors);
static void testnoisy(bool& idwerrors);

/*************************************************************************
Testing IDW interpolation
*************************************************************************/
bool testidw(bool silent)
{
    bool result;
    ap::real_2d_array xy;
    int i;
    int j;
    double vx;
    double vy;
    double vz;
    int d;
    int dtask;
    int nx;
    int n;
    int nq;
    int nw;
    int smalln;
    int largen;
    bool waserrors;
    bool idwerrors;

    idwerrors = false;
    smalln = 256;
    largen = 1024;
    nq = 10;
    nw = 18;
    
    //
    // Simple test:
    // * F = x^3 + sin(pi*y)*z^2 - (x+y)^2
    // * space is either R1=[-1,+1] (other dimensions are
    //   fixed at 0), R1^2 or R1^3.
    //* D = -1, 0, 1, 2
    //
    for(nx = 1; nx <= 2; nx++)
    {
        xy.setlength(largen, nx+1);
        for(i = 0; i <= largen-1; i++)
        {
            for(j = 0; j <= nx-1; j++)
            {
                xy(i,j) = 2*ap::randomreal()-1;
            }
            if( nx>=1 )
            {
                vx = xy(i,0);
            }
            else
            {
                vx = 0;
            }
            if( nx>=2 )
            {
                vy = xy(i,1);
            }
            else
            {
                vy = 0;
            }
            if( nx>=3 )
            {
                vz = xy(i,2);
            }
            else
            {
                vz = 0;
            }
            xy(i,nx) = vx*vx*vx+sin(ap::pi()*vy)*ap::sqr(vz)-ap::sqr(vx+vy);
        }
        for(d = -1; d <= 2; d++)
        {
            testxy(xy, largen, nx, d, nq, nw, idwerrors);
        }
    }
    
    //
    // Another simple test:
    // * five points in 2D - (0,0), (0,1), (1,0), (-1,0) (0,-1)
    // * F is random
    // * D = -1, 0, 1, 2
    //
    nx = 2;
    xy.setlength(5, nx+1);
    xy(0,0) = 0;
    xy(0,1) = 0;
    xy(0,2) = 2*ap::randomreal()-1;
    xy(1,0) = 1;
    xy(1,1) = 0;
    xy(1,2) = 2*ap::randomreal()-1;
    xy(2,0) = 0;
    xy(2,1) = 1;
    xy(2,2) = 2*ap::randomreal()-1;
    xy(3,0) = -1;
    xy(3,1) = 0;
    xy(3,2) = 2*ap::randomreal()-1;
    xy(4,0) = 0;
    xy(4,1) = -1;
    xy(4,2) = 2*ap::randomreal()-1;
    for(d = -1; d <= 2; d++)
    {
        testxy(xy, 5, nx, d, nq, nw, idwerrors);
    }
    
    //
    // Degree test.
    //
    // F is either:
    // * constant (DTask=0)
    // * linear (DTask=1)
    // * quadratic (DTask=2)
    //
    // Nodal functions are either
    // * constant (D=0)
    // * linear (D=1)
    // * quadratic (D=2)
    //
    // When DTask<=D, we can interpolate without errors.
    // When DTask>D, we MUST have errors.
    //
    for(nx = 1; nx <= 3; nx++)
    {
        for(d = 0; d <= 2; d++)
        {
            for(dtask = 0; dtask <= 2; dtask++)
            {
                testdegree(smalln, nx, d, dtask, idwerrors);
            }
        }
    }
    
    //
    // Noisy test
    //
    testnoisy(idwerrors);
    
    //
    // report
    //
    waserrors = idwerrors;
    if( !silent )
    {
        printf("TESTING INVERSE DISTANCE WEIGHTING\n");
        printf("* IDW:                                   ");
        if( !idwerrors )
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
Unsets 2D array.
*************************************************************************/
static void unset2d(ap::complex_2d_array& a)
{

    a.setbounds(0, 0, 0, 0);
    a(0,0) = 2*ap::randomreal()-1;
}


/*************************************************************************
Unsets 1D array.
*************************************************************************/
static void unset1d(ap::real_1d_array& a)
{

    a.setbounds(0, 0);
    a(0) = 2*ap::randomreal()-1;
}


/*************************************************************************
Testing IDW:
* generate model using N/NX/D/NQ/NW
* test basic properties
*************************************************************************/
static void testxy(const ap::real_2d_array& xy,
     const int& n,
     const int& nx,
     const int& d,
     const int& nq,
     const int& nw,
     bool& idwerrors)
{
    double threshold;
    double lipschitzstep;
    int i;
    int j;
    int i1;
    int i2;
    double v;
    double v1;
    double v2;
    double t;
    double l1;
    double l2;
    idwinterpolant z1;
    ap::real_1d_array x;

    threshold = 1000*ap::machineepsilon;
    lipschitzstep = 0.001;
    x.setlength(nx);
    
    //
    // build
    //
    idwbuildmodifiedshepard(xy, n, nx, d, nq, nw, z1);
    
    //
    // first, test interpolation properties at nodes
    //
    for(i = 0; i <= n-1; i++)
    {
        ap::vmove(&x(0), 1, &xy(i, 0), 1, ap::vlen(0,nx-1));
        idwerrors = idwerrors||ap::fp_neq(idwcalc(z1, x),xy(i,nx));
    }
    
    //
    // test Lipschitz continuity
    //
    i1 = ap::randominteger(n);
    do
    {
        i2 = ap::randominteger(n);
    }
    while(i2==i1);
    l1 = 0;
    t = 0;
    while(ap::fp_less(t,1))
    {
        v = 1-t;
        ap::vmove(&x(0), 1, &xy(i1, 0), 1, ap::vlen(0,nx-1), v);
        v = t;
        ap::vadd(&x(0), 1, &xy(i2, 0), 1, ap::vlen(0,nx-1), v);
        v1 = idwcalc(z1, x);
        v = 1-(t+lipschitzstep);
        ap::vmove(&x(0), 1, &xy(i1, 0), 1, ap::vlen(0,nx-1), v);
        v = t+lipschitzstep;
        ap::vadd(&x(0), 1, &xy(i2, 0), 1, ap::vlen(0,nx-1), v);
        v2 = idwcalc(z1, x);
        l1 = ap::maxreal(l1, fabs(v2-v1)/lipschitzstep);
        t = t+lipschitzstep;
    }
    l2 = 0;
    t = 0;
    while(ap::fp_less(t,1))
    {
        v = 1-t;
        ap::vmove(&x(0), 1, &xy(i1, 0), 1, ap::vlen(0,nx-1), v);
        v = t;
        ap::vadd(&x(0), 1, &xy(i2, 0), 1, ap::vlen(0,nx-1), v);
        v1 = idwcalc(z1, x);
        v = 1-(t+lipschitzstep/3);
        ap::vmove(&x(0), 1, &xy(i1, 0), 1, ap::vlen(0,nx-1), v);
        v = t+lipschitzstep/3;
        ap::vadd(&x(0), 1, &xy(i2, 0), 1, ap::vlen(0,nx-1), v);
        v2 = idwcalc(z1, x);
        l2 = ap::maxreal(l2, fabs(v2-v1)/(lipschitzstep/3));
        t = t+lipschitzstep/3;
    }
    idwerrors = idwerrors||ap::fp_greater(l2,2.0*l1);
}


/*************************************************************************
Testing IDW:
* generate model using R-based model
* test basic properties
*************************************************************************/
static void testrxy(const ap::real_2d_array& xy,
     const int& n,
     const int& nx,
     const double& r,
     bool& idwerrors)
{
    double threshold;
    double lipschitzstep;
    int i;
    int j;
    int i1;
    int i2;
    double v;
    double v1;
    double v2;
    double t;
    double l1;
    double l2;
    idwinterpolant z1;
    ap::real_1d_array x;

    threshold = 1000*ap::machineepsilon;
    lipschitzstep = 0.001;
    x.setlength(nx);
    
    //
    // build
    //
    idwbuildmodifiedshepardr(xy, n, nx, r, z1);
    
    //
    // first, test interpolation properties at nodes
    //
    for(i = 0; i <= n-1; i++)
    {
        ap::vmove(&x(0), 1, &xy(i, 0), 1, ap::vlen(0,nx-1));
        idwerrors = idwerrors||ap::fp_neq(idwcalc(z1, x),xy(i,nx));
    }
    
    //
    // test Lipschitz continuity
    //
    i1 = ap::randominteger(n);
    do
    {
        i2 = ap::randominteger(n);
    }
    while(i2==i1);
    l1 = 0;
    t = 0;
    while(ap::fp_less(t,1))
    {
        v = 1-t;
        ap::vmove(&x(0), 1, &xy(i1, 0), 1, ap::vlen(0,nx-1), v);
        v = t;
        ap::vadd(&x(0), 1, &xy(i2, 0), 1, ap::vlen(0,nx-1), v);
        v1 = idwcalc(z1, x);
        v = 1-(t+lipschitzstep);
        ap::vmove(&x(0), 1, &xy(i1, 0), 1, ap::vlen(0,nx-1), v);
        v = t+lipschitzstep;
        ap::vadd(&x(0), 1, &xy(i2, 0), 1, ap::vlen(0,nx-1), v);
        v2 = idwcalc(z1, x);
        l1 = ap::maxreal(l1, fabs(v2-v1)/lipschitzstep);
        t = t+lipschitzstep;
    }
    l2 = 0;
    t = 0;
    while(ap::fp_less(t,1))
    {
        v = 1-t;
        ap::vmove(&x(0), 1, &xy(i1, 0), 1, ap::vlen(0,nx-1), v);
        v = t;
        ap::vadd(&x(0), 1, &xy(i2, 0), 1, ap::vlen(0,nx-1), v);
        v1 = idwcalc(z1, x);
        v = 1-(t+lipschitzstep/3);
        ap::vmove(&x(0), 1, &xy(i1, 0), 1, ap::vlen(0,nx-1), v);
        v = t+lipschitzstep/3;
        ap::vadd(&x(0), 1, &xy(i2, 0), 1, ap::vlen(0,nx-1), v);
        v2 = idwcalc(z1, x);
        l2 = ap::maxreal(l2, fabs(v2-v1)/(lipschitzstep/3));
        t = t+lipschitzstep/3;
    }
    idwerrors = idwerrors||ap::fp_greater(l2,2.0*l1);
}


/*************************************************************************
Testing degree properties

F is either:
* constant (DTask=0)
* linear (DTask=1)
* quadratic (DTask=2)

Nodal functions are either
* constant (D=0)
* linear (D=1)
* quadratic (D=2)

When DTask<=D, we can interpolate without errors.
When DTask>D, we MUST have errors.
*************************************************************************/
static void testdegree(const int& n,
     const int& nx,
     const int& d,
     const int& dtask,
     bool& idwerrors)
{
    double threshold;
    int nq;
    int nw;
    int i;
    int j;
    double v;
    double c0;
    ap::real_1d_array c1;
    ap::real_2d_array c2;
    ap::real_1d_array x;
    ap::real_2d_array xy;
    idwinterpolant z1;
    double v1;
    double v2;

    threshold = 1.0E6*ap::machineepsilon;
    nq = 2*(nx*nx+nx+1);
    nw = 10;
    ap::ap_error::make_assertion(nq<=n, "TestDegree: internal error");
    
    //
    // prepare model
    //
    c0 = 2*ap::randomreal()-1;
    c1.setlength(nx);
    for(i = 0; i <= nx-1; i++)
    {
        c1(i) = 2*ap::randomreal()-1;
    }
    c2.setlength(nx, nx);
    for(i = 0; i <= nx-1; i++)
    {
        for(j = i+1; j <= nx-1; j++)
        {
            c2(i,j) = 2*ap::randomreal()-1;
            c2(j,i) = c2(i,j);
        }
        do
        {
            c2(i,i) = 2*ap::randomreal()-1;
        }
        while(ap::fp_less_eq(fabs(c2(i,i)),0.3));
    }
    
    //
    // prepare points
    //
    xy.setlength(n, nx+1);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= nx-1; j++)
        {
            xy(i,j) = 4*ap::randomreal()-2;
        }
        xy(i,nx) = c0;
        if( dtask>=1 )
        {
            v = ap::vdotproduct(&c1(0), 1, &xy(i, 0), 1, ap::vlen(0,nx-1));
            xy(i,nx) = xy(i,nx)+v;
        }
        if( dtask==2 )
        {
            for(j = 0; j <= nx-1; j++)
            {
                v = ap::vdotproduct(&c2(j, 0), 1, &xy(i, 0), 1, ap::vlen(0,nx-1));
                xy(i,nx) = xy(i,nx)+xy(i,j)*v;
            }
        }
    }
    
    //
    // build interpolant, calculate value at random point
    //
    idwbuildmodifiedshepard(xy, n, nx, d, nq, nw, z1);
    x.setlength(nx);
    for(i = 0; i <= nx-1; i++)
    {
        x(i) = 4*ap::randomreal()-2;
    }
    v1 = idwcalc(z1, x);
    
    //
    // calculate model value at the same point
    //
    v2 = c0;
    if( dtask>=1 )
    {
        v = ap::vdotproduct(&c1(0), 1, &x(0), 1, ap::vlen(0,nx-1));
        v2 = v2+v;
    }
    if( dtask==2 )
    {
        for(j = 0; j <= nx-1; j++)
        {
            v = ap::vdotproduct(&c2(j, 0), 1, &x(0), 1, ap::vlen(0,nx-1));
            v2 = v2+x(j)*v;
        }
    }
    
    //
    // Compare
    //
    if( dtask<=d )
    {
        idwerrors = idwerrors||ap::fp_greater(fabs(v2-v1),threshold);
    }
    else
    {
        idwerrors = idwerrors||ap::fp_less(fabs(v2-v1),threshold);
    }
}


/*************************************************************************
Noisy test:
 * F = x^2 + y^2 + z^2 + noise on [-1,+1]^3
 * space is either R1=[-1,+1] (other dimensions are
   fixed at 0), R1^2 or R1^3.
 * D = 1, 2
 * 4096 points is used for function generation,
   4096 points - for testing
 * RMS error of "noisy" model on test set must be
   lower than RMS error of interpolation model.
*************************************************************************/
static void testnoisy(bool& idwerrors)
{
    double noiselevel;
    int nq;
    int nw;
    int d;
    int nx;
    int ntrn;
    int ntst;
    int i;
    int j;
    double v;
    double t;
    double v1;
    double v2;
    double ve;
    ap::real_2d_array xy;
    ap::real_1d_array x;
    idwinterpolant z1;
    idwinterpolant z2;
    double rms1;
    double rms2;

    nq = 20;
    nw = 40;
    noiselevel = 0.2;
    ntrn = 256;
    ntst = 1024;
    for(d = 1; d <= 2; d++)
    {
        for(nx = 1; nx <= 2; nx++)
        {
            
            //
            // prepare dataset
            //
            xy.setlength(ntrn, nx+1);
            for(i = 0; i <= ntrn-1; i++)
            {
                v = noiselevel*(2*ap::randomreal()-1);
                for(j = 0; j <= nx-1; j++)
                {
                    t = 2*ap::randomreal()-1;
                    v = v+ap::sqr(t);
                    xy(i,j) = t;
                }
                xy(i,nx) = v;
            }
            
            //
            // build interpolants
            //
            idwbuildmodifiedshepard(xy, ntrn, nx, d, nq, nw, z1);
            idwbuildnoisy(xy, ntrn, nx, d, nq, nw, z2);
            
            //
            // calculate RMS errors
            //
            x.setlength(nx);
            rms1 = 0;
            rms2 = 0;
            for(i = 0; i <= ntst-1; i++)
            {
                ve = 0;
                for(j = 0; j <= nx-1; j++)
                {
                    t = 2*ap::randomreal()-1;
                    x(j) = t;
                    ve = ve+ap::sqr(t);
                }
                v1 = idwcalc(z1, x);
                v2 = idwcalc(z2, x);
                rms1 = rms1+ap::sqr(v1-ve);
                rms2 = rms2+ap::sqr(v2-ve);
            }
            idwerrors = idwerrors||ap::fp_greater(rms2,rms1);
        }
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testidwunit_test_silent()
{
    bool result;

    result = testidw(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testidwunit_test()
{
    bool result;

    result = testidw(false);
    return result;
}






 
#include <stdio.h>
#include "testpolintunit.h"

static double internalpolint(const ap::real_1d_array& x,
     ap::real_1d_array f,
     int n,
     double t);
static void brcunset(barycentricinterpolant& b);

/*************************************************************************
Unit test
*************************************************************************/
bool testpolint(bool silent)
{
    bool result;
    bool waserrors;
    bool interrors;
    bool fiterrors;
    double threshold;
    ap::real_1d_array x;
    ap::real_1d_array y;
    ap::real_1d_array w;
    ap::real_1d_array x2;
    ap::real_1d_array y2;
    ap::real_1d_array w2;
    ap::real_1d_array xfull;
    ap::real_1d_array yfull;
    double a;
    double b;
    double t;
    int i;
    int k;
    ap::real_1d_array xc;
    ap::real_1d_array yc;
    ap::integer_1d_array dc;
    int info;
    int info2;
    double v;
    double v0;
    double v1;
    double v2;
    double s;
    double xmin;
    double xmax;
    double refrms;
    double refavg;
    double refavgrel;
    double refmax;
    barycentricinterpolant p;
    barycentricinterpolant p1;
    barycentricinterpolant p2;
    polynomialfitreport rep;
    polynomialfitreport rep2;
    int n;
    int m;
    int maxn;
    int pass;
    int passcount;

    waserrors = false;
    interrors = false;
    fiterrors = false;
    maxn = 5;
    passcount = 20;
    threshold = 1.0E8*ap::machineepsilon;
    
    //
    // Test equidistant interpolation
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        for(n = 1; n <= maxn; n++)
        {
            
            //
            // prepare task:
            // * equidistant points
            // * random Y
            // * T in [A,B] or near (within 10% of its width)
            //
            do
            {
                a = 2*ap::randomreal()-1;
                b = 2*ap::randomreal()-1;
            }
            while(ap::fp_less_eq(fabs(a-b),0.2));
            t = a+(1.2*ap::randomreal()-0.1)*(b-a);
            taskgenint1dequidist(a, b, n, x, y);
            
            //
            // test "fast" equidistant interpolation (no barycentric model)
            //
            interrors = interrors||ap::fp_greater(fabs(polynomialcalceqdist(a, b, y, n, t)-internalpolint(x, y, n, t)),threshold);
            
            //
            // test "slow" equidistant interpolation (create barycentric model)
            //
            brcunset(p);
            polynomialbuild(x, y, n, p);
            interrors = interrors||ap::fp_greater(fabs(barycentriccalc(p, t)-internalpolint(x, y, n, t)),threshold);
            
            //
            // test "fast" interpolation (create "fast" barycentric model)
            //
            brcunset(p);
            polynomialbuildeqdist(a, b, y, n, p);
            interrors = interrors||ap::fp_greater(fabs(barycentriccalc(p, t)-internalpolint(x, y, n, t)),threshold);
        }
    }
    
    //
    // Test Chebyshev-1 interpolation
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        for(n = 1; n <= maxn; n++)
        {
            
            //
            // prepare task:
            // * equidistant points
            // * random Y
            // * T in [A,B] or near (within 10% of its width)
            //
            do
            {
                a = 2*ap::randomreal()-1;
                b = 2*ap::randomreal()-1;
            }
            while(ap::fp_less_eq(fabs(a-b),0.2));
            t = a+(1.2*ap::randomreal()-0.1)*(b-a);
            taskgenint1dcheb1(a, b, n, x, y);
            
            //
            // test "fast" interpolation (no barycentric model)
            //
            interrors = interrors||ap::fp_greater(fabs(polynomialcalccheb1(a, b, y, n, t)-internalpolint(x, y, n, t)),threshold);
            
            //
            // test "slow" interpolation (create barycentric model)
            //
            brcunset(p);
            polynomialbuild(x, y, n, p);
            interrors = interrors||ap::fp_greater(fabs(barycentriccalc(p, t)-internalpolint(x, y, n, t)),threshold);
            
            //
            // test "fast" interpolation (create "fast" barycentric model)
            //
            brcunset(p);
            polynomialbuildcheb1(a, b, y, n, p);
            interrors = interrors||ap::fp_greater(fabs(barycentriccalc(p, t)-internalpolint(x, y, n, t)),threshold);
        }
    }
    
    //
    // Test Chebyshev-2 interpolation
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        for(n = 1; n <= maxn; n++)
        {
            
            //
            // prepare task:
            // * equidistant points
            // * random Y
            // * T in [A,B] or near (within 10% of its width)
            //
            do
            {
                a = 2*ap::randomreal()-1;
                b = 2*ap::randomreal()-1;
            }
            while(ap::fp_less_eq(fabs(a-b),0.2));
            t = a+(1.2*ap::randomreal()-0.1)*(b-a);
            taskgenint1dcheb2(a, b, n, x, y);
            
            //
            // test "fast" interpolation (no barycentric model)
            //
            interrors = interrors||ap::fp_greater(fabs(polynomialcalccheb2(a, b, y, n, t)-internalpolint(x, y, n, t)),threshold);
            
            //
            // test "slow" interpolation (create barycentric model)
            //
            brcunset(p);
            polynomialbuild(x, y, n, p);
            interrors = interrors||ap::fp_greater(fabs(barycentriccalc(p, t)-internalpolint(x, y, n, t)),threshold);
            
            //
            // test "fast" interpolation (create "fast" barycentric model)
            //
            brcunset(p);
            polynomialbuildcheb2(a, b, y, n, p);
            interrors = interrors||ap::fp_greater(fabs(barycentriccalc(p, t)-internalpolint(x, y, n, t)),threshold);
        }
    }
    
    //
    // crash-test: ability to solve tasks which will overflow/underflow
    // weights with straightforward implementation
    //
    for(n = 1; n <= 20; n++)
    {
        a = -0.1*ap::maxrealnumber;
        b = +0.1*ap::maxrealnumber;
        taskgenint1dequidist(a, b, n, x, y);
        polynomialbuild(x, y, n, p);
        for(i = 0; i <= n-1; i++)
        {
            interrors = interrors||ap::fp_eq(p.w(i),0);
        }
    }
    
    //
    // Test rational fitting:
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        for(n = 1; n <= maxn; n++)
        {
            
            //
            // N=M+K fitting (i.e. interpolation)
            //
            for(k = 0; k <= n-1; k++)
            {
                taskgenint1d(double(-1), double(1), n, xfull, yfull);
                x.setlength(n-k);
                y.setlength(n-k);
                w.setlength(n-k);
                if( k>0 )
                {
                    xc.setlength(k);
                    yc.setlength(k);
                    dc.setlength(k);
                }
                for(i = 0; i <= n-k-1; i++)
                {
                    x(i) = xfull(i);
                    y(i) = yfull(i);
                    w(i) = 1+ap::randomreal();
                }
                for(i = 0; i <= k-1; i++)
                {
                    xc(i) = xfull(n-k+i);
                    yc(i) = yfull(n-k+i);
                    dc(i) = 0;
                }
                polynomialfitwc(x, y, w, n-k, xc, yc, dc, k, n, info, p1, rep);
                if( info<=0 )
                {
                    fiterrors = true;
                }
                else
                {
                    for(i = 0; i <= n-k-1; i++)
                    {
                        fiterrors = fiterrors||ap::fp_greater(fabs(barycentriccalc(p1, x(i))-y(i)),threshold);
                    }
                    for(i = 0; i <= k-1; i++)
                    {
                        fiterrors = fiterrors||ap::fp_greater(fabs(barycentriccalc(p1, xc(i))-yc(i)),threshold);
                    }
                }
            }
            
            //
            // Testing constraints on derivatives.
            // Special tasks which will always have solution:
            // 1. P(0)=YC[0]
            // 2. P(0)=YC[0], P'(0)=YC[1]
            //
            if( n>1 )
            {
                for(m = 3; m <= 5; m++)
                {
                    for(k = 1; k <= 2; k++)
                    {
                        taskgenint1d(double(-1), double(1), n, x, y);
                        w.setlength(n);
                        xc.setlength(2);
                        yc.setlength(2);
                        dc.setlength(2);
                        for(i = 0; i <= n-1; i++)
                        {
                            w(i) = 1+ap::randomreal();
                        }
                        xc(0) = 0;
                        yc(0) = 2*ap::randomreal()-1;
                        dc(0) = 0;
                        xc(1) = 0;
                        yc(1) = 2*ap::randomreal()-1;
                        dc(1) = 1;
                        polynomialfitwc(x, y, w, n, xc, yc, dc, k, m, info, p1, rep);
                        if( info<=0 )
                        {
                            fiterrors = true;
                        }
                        else
                        {
                            barycentricdiff1(p1, 0.0, v0, v1);
                            fiterrors = fiterrors||ap::fp_greater(fabs(v0-yc(0)),threshold);
                            if( k==2 )
                            {
                                fiterrors = fiterrors||ap::fp_greater(fabs(v1-yc(1)),threshold);
                            }
                        }
                    }
                }
            }
        }
    }
    for(m = 2; m <= 8; m++)
    {
        for(pass = 1; pass <= passcount; pass++)
        {
            
            //
            // General fitting
            //
            // interpolating function through M nodes should have
            // greater RMS error than fitting it through the same M nodes
            //
            n = 100;
            x2.setlength(n);
            y2.setlength(n);
            w2.setlength(n);
            xmin = 0;
            xmax = 2*ap::pi();
            for(i = 0; i <= n-1; i++)
            {
                x2(i) = 2*ap::pi()*ap::randomreal();
                y2(i) = sin(x2(i));
                w2(i) = 1;
            }
            x.setlength(m);
            y.setlength(m);
            for(i = 0; i <= m-1; i++)
            {
                x(i) = xmin+(xmax-xmin)*i/(m-1);
                y(i) = sin(x(i));
            }
            polynomialbuild(x, y, m, p1);
            polynomialfitwc(x2, y2, w2, n, xc, yc, dc, 0, m, info, p2, rep);
            if( info<=0 )
            {
                fiterrors = true;
            }
            else
            {
                
                //
                // calculate P1 (interpolant) RMS error, compare with P2 error
                //
                v1 = 0;
                v2 = 0;
                for(i = 0; i <= n-1; i++)
                {
                    v1 = v1+ap::sqr(barycentriccalc(p1, x2(i))-y2(i));
                    v2 = v2+ap::sqr(barycentriccalc(p2, x2(i))-y2(i));
                }
                v1 = sqrt(v1/n);
                v2 = sqrt(v2/n);
                fiterrors = fiterrors||ap::fp_greater(v2,v1);
                fiterrors = fiterrors||ap::fp_greater(fabs(v2-rep.rmserror),threshold);
            }
            
            //
            // compare weighted and non-weighted
            //
            n = 20;
            x.setlength(n);
            y.setlength(n);
            w.setlength(n);
            for(i = 0; i <= n-1; i++)
            {
                x(i) = 2*ap::randomreal()-1;
                y(i) = 2*ap::randomreal()-1;
                w(i) = 1;
            }
            polynomialfitwc(x, y, w, n, xc, yc, dc, 0, m, info, p1, rep);
            polynomialfit(x, y, n, m, info2, p2, rep2);
            if( info<=0||info2<=0 )
            {
                fiterrors = true;
            }
            else
            {
                
                //
                // calculate P1 (interpolant), compare with P2 error
                // compare RMS errors
                //
                t = 2*ap::randomreal()-1;
                v1 = barycentriccalc(p1, t);
                v2 = barycentriccalc(p2, t);
                fiterrors = fiterrors||ap::fp_neq(v2,v1);
                fiterrors = fiterrors||ap::fp_neq(rep.rmserror,rep2.rmserror);
                fiterrors = fiterrors||ap::fp_neq(rep.avgerror,rep2.avgerror);
                fiterrors = fiterrors||ap::fp_neq(rep.avgrelerror,rep2.avgrelerror);
                fiterrors = fiterrors||ap::fp_neq(rep.maxerror,rep2.maxerror);
            }
        }
    }
    for(pass = 1; pass <= passcount; pass++)
    {
        ap::ap_error::make_assertion(passcount>=2, "PassCount should be 2 or greater!");
        
        //
        // solve simple task (all X[] are the same, Y[] are specially
        // calculated to ensure simple form of all types of errors)
        // and check correctness of the errors calculated by subroutines
        //
        // First pass is done with zero Y[], other passes - with random Y[].
        // It should test both ability to correctly calculate errors and
        // ability to not fail while working with zeros :)
        //
        n = 4;
        if( pass==1 )
        {
            v1 = 0;
            v2 = 0;
            v = 0;
        }
        else
        {
            v1 = ap::randomreal();
            v2 = ap::randomreal();
            v = 1+ap::randomreal();
        }
        x.setlength(4);
        y.setlength(4);
        w.setlength(4);
        x(0) = 0;
        y(0) = v-v2;
        w(0) = 1;
        x(1) = 0;
        y(1) = v-v1;
        w(1) = 1;
        x(2) = 0;
        y(2) = v+v1;
        w(2) = 1;
        x(3) = 0;
        y(3) = v+v2;
        w(3) = 1;
        refrms = sqrt((ap::sqr(v1)+ap::sqr(v2))/2);
        refavg = (fabs(v1)+fabs(v2))/2;
        if( pass==1 )
        {
            refavgrel = 0;
        }
        else
        {
            refavgrel = 0.25*(fabs(v2)/fabs(v-v2)+fabs(v1)/fabs(v-v1)+fabs(v1)/fabs(v+v1)+fabs(v2)/fabs(v+v2));
        }
        refmax = ap::maxreal(v1, v2);
        
        //
        // Test errors correctness
        //
        polynomialfit(x, y, 4, 1, info, p, rep);
        if( info<=0 )
        {
            fiterrors = true;
        }
        else
        {
            s = barycentriccalc(p, double(0));
            fiterrors = fiterrors||ap::fp_greater(fabs(s-v),threshold);
            fiterrors = fiterrors||ap::fp_greater(fabs(rep.rmserror-refrms),threshold);
            fiterrors = fiterrors||ap::fp_greater(fabs(rep.avgerror-refavg),threshold);
            fiterrors = fiterrors||ap::fp_greater(fabs(rep.avgrelerror-refavgrel),threshold);
            fiterrors = fiterrors||ap::fp_greater(fabs(rep.maxerror-refmax),threshold);
        }
    }
    
    //
    // report
    //
    waserrors = interrors||fiterrors;
    if( !silent )
    {
        printf("TESTING POLYNOMIAL INTERPOLATION AND FITTING\n");
        
        //
        // Normal tests
        //
        printf("INTERPOLATION TEST:                      ");
        if( interrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("FITTING TEST:                            ");
        if( fiterrors )
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
        printf("\n\n");
    }
    
    //
    // end
    //
    result = !waserrors;
    return result;
}


static double internalpolint(const ap::real_1d_array& x,
     ap::real_1d_array f,
     int n,
     double t)
{
    double result;
    int i;
    int j;

    n = n-1;
    for(j = 0; j <= n-1; j++)
    {
        for(i = j+1; i <= n; i++)
        {
            f(i) = ((t-x(j))*f(i)-(t-x(i))*f(j))/(x(i)-x(j));
        }
    }
    result = f(n);
    return result;
}


static void brcunset(barycentricinterpolant& b)
{
    ap::real_1d_array x;
    ap::real_1d_array y;
    ap::real_1d_array w;

    x.setlength(1);
    y.setlength(1);
    w.setlength(1);
    x(0) = 0;
    y(0) = 0;
    w(0) = 1;
    barycentricbuildxyw(x, y, w, 1, b);
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testpolintunit_test_silent()
{
    bool result;

    result = testpolint(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testpolintunit_test()
{
    bool result;

    result = testpolint(false);
    return result;
}





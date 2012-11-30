
 
#include <stdio.h>
#include "testratinterpolation.h"

static void poldiff2(const ap::real_1d_array& x,
     ap::real_1d_array f,
     int n,
     double t,
     double& p,
     double& dp,
     double& d2p);
static void brcunset(barycentricinterpolant& b);
static bool is1dsolution(int n,
     const ap::real_1d_array& y,
     const ap::real_1d_array& w,
     double c);

bool testri(bool silent)
{
    bool result;
    bool waserrors;
    bool bcerrors;
    bool nperrors;
    bool fiterrors;
    double threshold;
    double lipschitztol;
    int maxn;
    int passcount;
    barycentricinterpolant b1;
    barycentricinterpolant b2;
    ap::real_1d_array x;
    ap::real_1d_array x2;
    ap::real_1d_array y;
    ap::real_1d_array y2;
    ap::real_1d_array w;
    ap::real_1d_array w2;
    ap::real_1d_array xc;
    ap::real_1d_array yc;
    ap::integer_1d_array dc;
    double h;
    double s1;
    double s2;
    bool bsame;
    int n;
    int m;
    int n2;
    int i;
    int j;
    int k;
    int d;
    int pass;
    double err;
    double maxerr;
    double t;
    double a;
    double b;
    double s;
    double v;
    double v0;
    double v1;
    double v2;
    double v3;
    double d0;
    double d1;
    double d2;
    int info;
    int info2;
    double xmin;
    double xmax;
    double refrms;
    double refavg;
    double refavgrel;
    double refmax;
    ap::real_1d_array ra;
    ap::real_1d_array ra2;
    int ralen;
    barycentricfitreport rep;
    barycentricfitreport rep2;
    barycentricinterpolant b3;
    barycentricinterpolant b4;

    nperrors = false;
    bcerrors = false;
    fiterrors = false;
    waserrors = false;
    
    //
    // PassCount        number of repeated passes
    // Threshold        error tolerance
    // LipschitzTol     Lipschitz constant increase allowed
    //                  when calculating constant on a twice denser grid
    //
    passcount = 5;
    maxn = 15;
    threshold = 1000000*ap::machineepsilon;
    lipschitztol = 1.3;
    
    //
    // Basic barycentric functions
    //
    for(n = 1; n <= 10; n++)
    {
        
        //
        // randomized tests
        //
        for(pass = 1; pass <= passcount; pass++)
        {
            
            //
            // generate weights from polynomial interpolation
            //
            v0 = 1+0.4*ap::randomreal()-0.2;
            v1 = 2*ap::randomreal()-1;
            v2 = 2*ap::randomreal()-1;
            v3 = 2*ap::randomreal()-1;
            x.setlength(n);
            y.setlength(n);
            w.setlength(n);
            for(i = 0; i <= n-1; i++)
            {
                if( n==1 )
                {
                    x(i) = 0;
                }
                else
                {
                    x(i) = v0*cos(i*ap::pi()/(n-1));
                }
                y(i) = sin(v1*x(i))+cos(v2*x(i))+exp(v3*x(i));
            }
            for(j = 0; j <= n-1; j++)
            {
                w(j) = 1;
                for(k = 0; k <= n-1; k++)
                {
                    if( k!=j )
                    {
                        w(j) = w(j)/(x(j)-x(k));
                    }
                }
            }
            barycentricbuildxyw(x, y, w, n, b1);
            
            //
            // unpack, then pack again and compare
            //
            brcunset(b2);
            barycentricunpack(b1, n2, x2, y2, w2);
            bcerrors = bcerrors||n2!=n;
            barycentricbuildxyw(x2, y2, w2, n2, b2);
            t = 2*ap::randomreal()-1;
            bcerrors = bcerrors||ap::fp_greater(fabs(barycentriccalc(b1, t)-barycentriccalc(b2, t)),threshold);
            
            //
            // serialize, unserialize, compare
            //
            brcunset(b2);
            barycentricserialize(b1, ra, ralen);
            ra2.setlength(ralen);
            ap::vmove(&ra2(0), 1, &ra(0), 1, ap::vlen(0,ralen-1));
            barycentricunserialize(ra2, b2);
            t = 2*ap::randomreal()-1;
            bcerrors = bcerrors||ap::fp_greater(fabs(barycentriccalc(b1, t)-barycentriccalc(b2, t)),threshold);
            
            //
            // copy, compare
            //
            brcunset(b2);
            barycentriccopy(b1, b2);
            t = 2*ap::randomreal()-1;
            bcerrors = bcerrors||ap::fp_greater(fabs(barycentriccalc(b1, t)-barycentriccalc(b2, t)),threshold);
            
            //
            // test interpolation properties
            //
            for(i = 0; i <= n-1; i++)
            {
                
                //
                // test interpolation at nodes
                //
                bcerrors = bcerrors||ap::fp_greater(fabs(barycentriccalc(b1, x(i))-y(i)),threshold*fabs(y(i)));
                
                //
                // compare with polynomial interpolation
                //
                t = 2*ap::randomreal()-1;
                poldiff2(x, y, n, t, v0, v1, v2);
                bcerrors = bcerrors||ap::fp_greater(fabs(barycentriccalc(b1, t)-v0),threshold*ap::maxreal(fabs(v0), double(1)));
                
                //
                // test continuity between nodes
                // calculate Lipschitz constant on two grids -
                // dense and even more dense. If Lipschitz constant
                // on a denser grid is significantly increased,
                // continuity test is failed
                //
                t = 3.0;
                k = 100;
                s1 = 0;
                for(j = 0; j <= k-1; j++)
                {
                    v1 = x(i)+(t-x(i))*j/k;
                    v2 = x(i)+(t-x(i))*(j+1)/k;
                    s1 = ap::maxreal(s1, fabs(barycentriccalc(b1, v2)-barycentriccalc(b1, v1))/fabs(v2-v1));
                }
                k = 2*k;
                s2 = 0;
                for(j = 0; j <= k-1; j++)
                {
                    v1 = x(i)+(t-x(i))*j/k;
                    v2 = x(i)+(t-x(i))*(j+1)/k;
                    s2 = ap::maxreal(s2, fabs(barycentriccalc(b1, v2)-barycentriccalc(b1, v1))/fabs(v2-v1));
                }
                bcerrors = bcerrors||ap::fp_greater(s2,lipschitztol*s1)&&ap::fp_greater(s1,threshold*k);
            }
            
            //
            // test differentiation properties
            //
            for(i = 0; i <= n-1; i++)
            {
                t = 2*ap::randomreal()-1;
                poldiff2(x, y, n, t, v0, v1, v2);
                d0 = 0;
                d1 = 0;
                d2 = 0;
                barycentricdiff1(b1, t, d0, d1);
                bcerrors = bcerrors||ap::fp_greater(fabs(v0-d0),threshold*ap::maxreal(fabs(v0), double(1)));
                bcerrors = bcerrors||ap::fp_greater(fabs(v1-d1),threshold*ap::maxreal(fabs(v1), double(1)));
                d0 = 0;
                d1 = 0;
                d2 = 0;
                barycentricdiff2(b1, t, d0, d1, d2);
                bcerrors = bcerrors||ap::fp_greater(fabs(v0-d0),threshold*ap::maxreal(fabs(v0), double(1)));
                bcerrors = bcerrors||ap::fp_greater(fabs(v1-d1),threshold*ap::maxreal(fabs(v1), double(1)));
                bcerrors = bcerrors||ap::fp_greater(fabs(v2-d2),sqrt(threshold)*ap::maxreal(fabs(v2), double(1)));
            }
            
            //
            // test linear translation
            //
            t = 2*ap::randomreal()-1;
            a = 2*ap::randomreal()-1;
            b = 2*ap::randomreal()-1;
            brcunset(b2);
            barycentriccopy(b1, b2);
            barycentriclintransx(b2, a, b);
            bcerrors = bcerrors||ap::fp_greater(fabs(barycentriccalc(b1, a*t+b)-barycentriccalc(b2, t)),threshold);
            a = 0;
            b = 2*ap::randomreal()-1;
            brcunset(b2);
            barycentriccopy(b1, b2);
            barycentriclintransx(b2, a, b);
            bcerrors = bcerrors||ap::fp_greater(fabs(barycentriccalc(b1, a*t+b)-barycentriccalc(b2, t)),threshold);
            a = 2*ap::randomreal()-1;
            b = 2*ap::randomreal()-1;
            brcunset(b2);
            barycentriccopy(b1, b2);
            barycentriclintransy(b2, a, b);
            bcerrors = bcerrors||ap::fp_greater(fabs(a*barycentriccalc(b1, t)+b-barycentriccalc(b2, t)),threshold);
        }
    }
    for(pass = 0; pass <= 3; pass++)
    {
        
        //
        // Crash-test: small numbers, large numbers
        //
        x.setlength(4);
        y.setlength(4);
        w.setlength(4);
        h = 1;
        if( pass%2==0 )
        {
            h = 100*ap::minrealnumber;
        }
        if( pass%2==1 )
        {
            h = 0.01*ap::maxrealnumber;
        }
        x(0) = 0*h;
        x(1) = 1*h;
        x(2) = 2*h;
        x(3) = 3*h;
        y(0) = 0*h;
        y(1) = 1*h;
        y(2) = 2*h;
        y(3) = 3*h;
        w(0) = -1/(x(1)-x(0));
        w(1) = +1*(1/(x(1)-x(0))+1/(x(2)-x(1)));
        w(2) = -1*(1/(x(2)-x(1))+1/(x(3)-x(2)));
        w(3) = +1/(x(3)-x(2));
        if( pass/2==0 )
        {
            v0 = 0;
        }
        if( pass/2==1 )
        {
            v0 = 0.6*h;
        }
        barycentricbuildxyw(x, y, w, 4, b1);
        t = barycentriccalc(b1, v0);
        d0 = 0;
        d1 = 0;
        d2 = 0;
        barycentricdiff1(b1, v0, d0, d1);
        bcerrors = bcerrors||ap::fp_greater(fabs(t-v0),threshold*v0);
        bcerrors = bcerrors||ap::fp_greater(fabs(d0-v0),threshold*v0);
        bcerrors = bcerrors||ap::fp_greater(fabs(d1-1),1000*threshold);
    }
    
    //
    // crash test: large abscissas, small argument
    //
    // test for errors in D0 is not very strict
    // because renormalization used in Diff1()
    // destroys part of precision.
    //
    x.setlength(4);
    y.setlength(4);
    w.setlength(4);
    h = 0.01*ap::maxrealnumber;
    x(0) = 0*h;
    x(1) = 1*h;
    x(2) = 2*h;
    x(3) = 3*h;
    y(0) = 0*h;
    y(1) = 1*h;
    y(2) = 2*h;
    y(3) = 3*h;
    w(0) = -1/(x(1)-x(0));
    w(1) = +1*(1/(x(1)-x(0))+1/(x(2)-x(1)));
    w(2) = -1*(1/(x(2)-x(1))+1/(x(3)-x(2)));
    w(3) = +1/(x(3)-x(2));
    v0 = 100*ap::minrealnumber;
    barycentricbuildxyw(x, y, w, 4, b1);
    t = barycentriccalc(b1, v0);
    d0 = 0;
    d1 = 0;
    d2 = 0;
    barycentricdiff1(b1, v0, d0, d1);
    bcerrors = bcerrors||ap::fp_greater(fabs(t),v0*(1+threshold));
    bcerrors = bcerrors||ap::fp_greater(fabs(d0),v0*(1+threshold));
    bcerrors = bcerrors||ap::fp_greater(fabs(d1-1),1000*threshold);
    
    //
    // crash test: test safe barycentric formula
    //
    x.setlength(4);
    y.setlength(4);
    w.setlength(4);
    h = 2*ap::minrealnumber;
    x(0) = 0*h;
    x(1) = 1*h;
    x(2) = 2*h;
    x(3) = 3*h;
    y(0) = 0*h;
    y(1) = 1*h;
    y(2) = 2*h;
    y(3) = 3*h;
    w(0) = -1/(x(1)-x(0));
    w(1) = +1*(1/(x(1)-x(0))+1/(x(2)-x(1)));
    w(2) = -1*(1/(x(2)-x(1))+1/(x(3)-x(2)));
    w(3) = +1/(x(3)-x(2));
    v0 = ap::minrealnumber;
    barycentricbuildxyw(x, y, w, 4, b1);
    t = barycentriccalc(b1, v0);
    bcerrors = bcerrors||ap::fp_greater(fabs(t-v0)/v0,threshold);
    
    //
    // Testing "No Poles" interpolation
    //
    maxerr = 0;
    for(pass = 1; pass <= passcount-1; pass++)
    {
        x.setlength(1);
        y.setlength(1);
        x(0) = 2*ap::randomreal()-1;
        y(0) = 2*ap::randomreal()-1;
        barycentricbuildfloaterhormann(x, y, 1, 1, b1);
        maxerr = ap::maxreal(maxerr, fabs(barycentriccalc(b1, 2*ap::randomreal()-1)-y(0)));
    }
    for(n = 2; n <= 10; n++)
    {
        
        //
        // compare interpolant built by subroutine
        // with interpolant built by hands
        //
        x.setlength(n);
        y.setlength(n);
        w.setlength(n);
        w2.setlength(n);
        
        //
        // D=1, non-equidistant nodes
        //
        for(pass = 1; pass <= passcount; pass++)
        {
            
            //
            // Initialize X, Y, W
            //
            a = -1-1*ap::randomreal();
            b = +1+1*ap::randomreal();
            for(i = 0; i <= n-1; i++)
            {
                x(i) = atan((b-a)*i/(n-1)+a);
            }
            for(i = 0; i <= n-1; i++)
            {
                y(i) = 2*ap::randomreal()-1;
            }
            w(0) = -1/(x(1)-x(0));
            s = 1;
            for(i = 1; i <= n-2; i++)
            {
                w(i) = s*(1/(x(i)-x(i-1))+1/(x(i+1)-x(i)));
                s = -s;
            }
            w(n-1) = s/(x(n-1)-x(n-2));
            for(i = 0; i <= n-1; i++)
            {
                k = ap::randominteger(n);
                if( k!=i )
                {
                    t = x(i);
                    x(i) = x(k);
                    x(k) = t;
                    t = y(i);
                    y(i) = y(k);
                    y(k) = t;
                    t = w(i);
                    w(i) = w(k);
                    w(k) = t;
                }
            }
            
            //
            // Build and test
            //
            barycentricbuildfloaterhormann(x, y, n, 1, b1);
            barycentricbuildxyw(x, y, w, n, b2);
            for(i = 1; i <= 2*n; i++)
            {
                t = a+(b-a)*ap::randomreal();
                maxerr = ap::maxreal(maxerr, fabs(barycentriccalc(b1, t)-barycentriccalc(b2, t)));
            }
        }
        
        //
        // D = 0, 1, 2. Equidistant nodes.
        //
        for(d = 0; d <= 2; d++)
        {
            for(pass = 1; pass <= passcount; pass++)
            {
                
                //
                // Skip incorrect (N,D) pairs
                //
                if( n<2*d )
                {
                    continue;
                }
                
                //
                // Initialize X, Y, W
                //
                a = -1-1*ap::randomreal();
                b = +1+1*ap::randomreal();
                for(i = 0; i <= n-1; i++)
                {
                    x(i) = (b-a)*i/(n-1)+a;
                }
                for(i = 0; i <= n-1; i++)
                {
                    y(i) = 2*ap::randomreal()-1;
                }
                s = 1;
                if( d==0 )
                {
                    for(i = 0; i <= n-1; i++)
                    {
                        w(i) = s;
                        s = -s;
                    }
                }
                if( d==1 )
                {
                    w(0) = -s;
                    for(i = 1; i <= n-2; i++)
                    {
                        w(i) = 2*s;
                        s = -s;
                    }
                    w(n-1) = s;
                }
                if( d==2 )
                {
                    w(0) = s;
                    w(1) = -3*s;
                    for(i = 2; i <= n-3; i++)
                    {
                        w(i) = 4*s;
                        s = -s;
                    }
                    w(n-2) = 3*s;
                    w(n-1) = -s;
                }
                
                //
                // Mix
                //
                for(i = 0; i <= n-1; i++)
                {
                    k = ap::randominteger(n);
                    if( k!=i )
                    {
                        t = x(i);
                        x(i) = x(k);
                        x(k) = t;
                        t = y(i);
                        y(i) = y(k);
                        y(k) = t;
                        t = w(i);
                        w(i) = w(k);
                        w(k) = t;
                    }
                }
                
                //
                // Build and test
                //
                barycentricbuildfloaterhormann(x, y, n, d, b1);
                barycentricbuildxyw(x, y, w, n, b2);
                for(i = 1; i <= 2*n; i++)
                {
                    t = a+(b-a)*ap::randomreal();
                    maxerr = ap::maxreal(maxerr, fabs(barycentriccalc(b1, t)-barycentriccalc(b2, t)));
                }
            }
        }
    }
    if( ap::fp_greater(maxerr,threshold) )
    {
        nperrors = true;
    }
    
    //
    // Test rational fitting:
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        for(n = 2; n <= maxn; n++)
        {
            
            //
            // N=M+K fitting (i.e. interpolation)
            //
            for(k = 0; k <= n-1; k++)
            {
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
                    x(i) = double(i)/double(n-1);
                    y(i) = 2*ap::randomreal()-1;
                    w(i) = 1+ap::randomreal();
                }
                for(i = 0; i <= k-1; i++)
                {
                    xc(i) = double(n-k+i)/double(n-1);
                    yc(i) = 2*ap::randomreal()-1;
                    dc(i) = 0;
                }
                barycentricfitfloaterhormannwc(x, y, w, n-k, xc, yc, dc, k, n, info, b1, rep);
                if( info<=0 )
                {
                    fiterrors = true;
                }
                else
                {
                    for(i = 0; i <= n-k-1; i++)
                    {
                        fiterrors = fiterrors||ap::fp_greater(fabs(barycentriccalc(b1, x(i))-y(i)),threshold);
                    }
                    for(i = 0; i <= k-1; i++)
                    {
                        fiterrors = fiterrors||ap::fp_greater(fabs(barycentriccalc(b1, xc(i))-yc(i)),threshold);
                    }
                }
            }
            
            //
            // Testing constraints on derivatives:
            // * several M's are tried
            // * several K's are tried - 1, 2.
            // * constraints at the ends of the interval
            //
            for(m = 3; m <= 5; m++)
            {
                for(k = 1; k <= 2; k++)
                {
                    x.setlength(n);
                    y.setlength(n);
                    w.setlength(n);
                    xc.setlength(2);
                    yc.setlength(2);
                    dc.setlength(2);
                    for(i = 0; i <= n-1; i++)
                    {
                        x(i) = 2*ap::randomreal()-1;
                        y(i) = 2*ap::randomreal()-1;
                        w(i) = 1+ap::randomreal();
                    }
                    xc(0) = -1;
                    yc(0) = 2*ap::randomreal()-1;
                    dc(0) = 0;
                    xc(1) = +1;
                    yc(1) = 2*ap::randomreal()-1;
                    dc(1) = 0;
                    barycentricfitfloaterhormannwc(x, y, w, n, xc, yc, dc, k, m, info, b1, rep);
                    if( info<=0 )
                    {
                        fiterrors = true;
                    }
                    else
                    {
                        for(i = 0; i <= k-1; i++)
                        {
                            barycentricdiff1(b1, xc(i), v0, v1);
                            fiterrors = fiterrors||ap::fp_greater(fabs(v0-yc(i)),threshold);
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
            xmin = ap::maxrealnumber;
            xmax = -ap::maxrealnumber;
            for(i = 0; i <= n-1; i++)
            {
                x2(i) = 2*ap::pi()*ap::randomreal();
                y2(i) = sin(x2(i));
                w2(i) = 1;
                xmin = ap::minreal(xmin, x2(i));
                xmax = ap::maxreal(xmax, x2(i));
            }
            x.setlength(m);
            y.setlength(m);
            for(i = 0; i <= m-1; i++)
            {
                x(i) = xmin+(xmax-xmin)*i/(m-1);
                y(i) = sin(x(i));
            }
            barycentricbuildfloaterhormann(x, y, m, 3, b1);
            barycentricfitfloaterhormannwc(x2, y2, w2, n, xc, yc, dc, 0, m, info, b2, rep);
            if( info<=0 )
            {
                fiterrors = true;
            }
            else
            {
                
                //
                // calculate B1 (interpolant) RMS error, compare with B2 error
                //
                v1 = 0;
                v2 = 0;
                for(i = 0; i <= n-1; i++)
                {
                    v1 = v1+ap::sqr(barycentriccalc(b1, x2(i))-y2(i));
                    v2 = v2+ap::sqr(barycentriccalc(b2, x2(i))-y2(i));
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
            barycentricfitfloaterhormannwc(x, y, w, n, xc, yc, dc, 0, m, info, b1, rep);
            barycentricfitfloaterhormann(x, y, n, m, info2, b2, rep2);
            if( info<=0||info2<=0 )
            {
                fiterrors = true;
            }
            else
            {
                
                //
                // calculate B1 (interpolant), compare with B2
                // compare RMS errors
                //
                t = 2*ap::randomreal()-1;
                v1 = barycentriccalc(b1, t);
                v2 = barycentriccalc(b2, t);
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
        barycentricfitfloaterhormann(x, y, 4, 2, info, b1, rep);
        if( info<=0 )
        {
            fiterrors = true;
        }
        else
        {
            s = barycentriccalc(b1, double(0));
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
    waserrors = bcerrors||nperrors||fiterrors;
    if( !silent )
    {
        printf("TESTING RATIONAL INTERPOLATION\n");
        printf("BASIC BARYCENTRIC FUNCTIONS:             ");
        if( bcerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("FLOATER-HORMANN:                         ");
        if( nperrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("RATIONAL FITTING:                        ");
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


static void poldiff2(const ap::real_1d_array& x,
     ap::real_1d_array f,
     int n,
     double t,
     double& p,
     double& dp,
     double& d2p)
{
    int m;
    int i;
    ap::real_1d_array df;
    ap::real_1d_array d2f;

    n = n-1;
    df.setbounds(0, n);
    d2f.setbounds(0, n);
    for(i = 0; i <= n; i++)
    {
        d2f(i) = 0;
        df(i) = 0;
    }
    for(m = 1; m <= n; m++)
    {
        for(i = 0; i <= n-m; i++)
        {
            d2f(i) = ((t-x(i+m))*d2f(i)+(x(i)-t)*d2f(i+1)+2*df(i)-2*df(i+1))/(x(i)-x(i+m));
            df(i) = ((t-x(i+m))*df(i)+f(i)+(x(i)-t)*df(i+1)-f(i+1))/(x(i)-x(i+m));
            f(i) = ((t-x(i+m))*f(i)+(x(i)-t)*f(i+1))/(x(i)-x(i+m));
        }
    }
    p = f(0);
    dp = df(0);
    d2p = d2f(0);
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
Tests whether constant C is solution of 1D LLS problem
*************************************************************************/
static bool is1dsolution(int n,
     const ap::real_1d_array& y,
     const ap::real_1d_array& w,
     double c)
{
    bool result;
    int i;
    double v;
    double s1;
    double s2;
    double s3;
    double delta;

    delta = 0.001;
    
    //
    // Test result
    //
    s1 = 0;
    for(i = 0; i <= n-1; i++)
    {
        s1 = s1+ap::sqr(w(i)*(c-y(i)));
    }
    s2 = 0;
    s3 = 0;
    for(i = 0; i <= n-1; i++)
    {
        s2 = s2+ap::sqr(w(i)*(c+delta-y(i)));
        s3 = s3+ap::sqr(w(i)*(c-delta-y(i)));
    }
    result = ap::fp_greater_eq(s2,s1)&&ap::fp_greater_eq(s3,s1);
    return result;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testratinterpolation_test_silent()
{
    bool result;

    result = testri(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testratinterpolation_test()
{
    bool result;

    result = testri(false);
    return result;
}





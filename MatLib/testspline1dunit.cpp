
 
#include <stdio.h>
#include "testspline1dunit.h"

static void lconst(double a,
     double b,
     const spline1dinterpolant& c,
     double lstep,
     double& l0,
     double& l1,
     double& l2);
static bool testunpack(const spline1dinterpolant& c,
     const ap::real_1d_array& x);
static void unsetspline1d(spline1dinterpolant& c);
static void unset1d(ap::real_1d_array& x);
static bool is1dsolution(int n,
     const ap::real_1d_array& y,
     const ap::real_1d_array& w,
     double c);

bool testsplineinterpolation(bool silent)
{
    bool result;
    bool waserrors;
    bool cserrors;
    bool hserrors;
    bool aserrors;
    bool lserrors;
    bool dserrors;
    bool uperrors;
    bool cperrors;
    bool lterrors;
    bool ierrors;
    bool fiterrors;
    double nonstrictthreshold;
    double threshold;
    int passcount;
    double lstep;
    double h;
    int maxn;
    int n;
    int m;
    int i;
    int k;
    int pass;
    int bltype;
    int brtype;
    int stype;
    ap::real_1d_array x;
    ap::real_1d_array y;
    ap::real_1d_array w;
    ap::real_1d_array w2;
    ap::real_1d_array y2;
    ap::real_1d_array d;
    ap::real_1d_array xc;
    ap::real_1d_array yc;
    ap::integer_1d_array dc;
    spline1dinterpolant c;
    spline1dinterpolant c2;
    int info;
    int info1;
    int info2;
    double a;
    double b;
    double bl;
    double br;
    double t;
    double sa;
    double sb;
    double v;
    double v1;
    double v2;
    double l10;
    double l11;
    double l12;
    double l20;
    double l21;
    double l22;
    double p0;
    double p1;
    double p2;
    double s;
    double ds;
    double d2s;
    double s2;
    double ds2;
    double d2s2;
    double vl;
    double vm;
    double vr;
    double err;
    spline1dfitreport rep;
    spline1dfitreport rep2;
    double refrms;
    double refavg;
    double refavgrel;
    double refmax;
    ap::real_1d_array ra;
    ap::real_1d_array ra2;
    int ralen;

    waserrors = false;
    passcount = 20;
    lstep = 0.005;
    h = 0.00001;
    maxn = 10;
    threshold = 10000*ap::machineepsilon;
    nonstrictthreshold = 0.00001;
    lserrors = false;
    cserrors = false;
    hserrors = false;
    aserrors = false;
    dserrors = false;
    cperrors = false;
    uperrors = false;
    lterrors = false;
    ierrors = false;
    fiterrors = false;
    
    //
    // General test: linear, cubic, Hermite, Akima
    //
    for(n = 2; n <= maxn; n++)
    {
        x.setbounds(0, n-1);
        y.setbounds(0, n-1);
        d.setbounds(0, n-1);
        for(pass = 1; pass <= passcount; pass++)
        {
            
            //
            // Prepare task
            //
            a = -1-ap::randomreal();
            b = +1+ap::randomreal();
            bl = 2*ap::randomreal()-1;
            br = 2*ap::randomreal()-1;
            for(i = 0; i <= n-1; i++)
            {
                x(i) = 0.5*(b+a)+0.5*(b-a)*cos(ap::pi()*(2*i+1)/(2*n));
                if( i==0 )
                {
                    x(i) = a;
                }
                if( i==n-1 )
                {
                    x(i) = b;
                }
                y(i) = cos(1.3*ap::pi()*x(i)+0.4);
                d(i) = -1.3*ap::pi()*sin(1.3*ap::pi()*x(i)+0.4);
            }
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
                    t = d(i);
                    d(i) = d(k);
                    d(k) = t;
                }
            }
            
            //
            // Build linear spline
            // Test for general interpolation scheme properties:
            // * values at nodes
            // * continuous function
            // Test for specific properties is implemented below.
            //
            spline1dbuildlinear(x, y, n, c);
            err = 0;
            for(i = 0; i <= n-1; i++)
            {
                err = ap::maxreal(err, fabs(y(i)-spline1dcalc(c, x(i))));
            }
            lserrors = lserrors||ap::fp_greater(err,threshold);
            lconst(a, b, c, lstep, l10, l11, l12);
            lconst(a, b, c, lstep/3, l20, l21, l22);
            lserrors = lserrors||ap::fp_greater(l20/l10,1.2);
            
            //
            // Build cubic spline.
            // Test for interpolation scheme properties:
            // * values at nodes
            // * boundary conditions
            // * continuous function
            // * continuous first derivative
            // * continuous second derivative
            //
            for(bltype = 0; bltype <= 2; bltype++)
            {
                for(brtype = 0; brtype <= 2; brtype++)
                {
                    spline1dbuildcubic(x, y, n, bltype, bl, brtype, br, c);
                    err = 0;
                    for(i = 0; i <= n-1; i++)
                    {
                        err = ap::maxreal(err, fabs(y(i)-spline1dcalc(c, x(i))));
                    }
                    cserrors = cserrors||ap::fp_greater(err,threshold);
                    err = 0;
                    if( bltype==0 )
                    {
                        spline1ddiff(c, a-h, s, ds, d2s);
                        spline1ddiff(c, a+h, s2, ds2, d2s2);
                        t = (d2s2-d2s)/(2*h);
                        err = ap::maxreal(err, fabs(t));
                    }
                    if( bltype==1 )
                    {
                        t = (spline1dcalc(c, a+h)-spline1dcalc(c, a-h))/(2*h);
                        err = ap::maxreal(err, fabs(bl-t));
                    }
                    if( bltype==2 )
                    {
                        t = (spline1dcalc(c, a+h)-2*spline1dcalc(c, a)+spline1dcalc(c, a-h))/ap::sqr(h);
                        err = ap::maxreal(err, fabs(bl-t));
                    }
                    if( brtype==0 )
                    {
                        spline1ddiff(c, b-h, s, ds, d2s);
                        spline1ddiff(c, b+h, s2, ds2, d2s2);
                        t = (d2s2-d2s)/(2*h);
                        err = ap::maxreal(err, fabs(t));
                    }
                    if( brtype==1 )
                    {
                        t = (spline1dcalc(c, b+h)-spline1dcalc(c, b-h))/(2*h);
                        err = ap::maxreal(err, fabs(br-t));
                    }
                    if( brtype==2 )
                    {
                        t = (spline1dcalc(c, b+h)-2*spline1dcalc(c, b)+spline1dcalc(c, b-h))/ap::sqr(h);
                        err = ap::maxreal(err, fabs(br-t));
                    }
                    cserrors = cserrors||ap::fp_greater(err,1.0E-3);
                    lconst(a, b, c, lstep, l10, l11, l12);
                    lconst(a, b, c, lstep/3, l20, l21, l22);
                    cserrors = cserrors||ap::fp_greater(l20/l10,1.2)&&ap::fp_greater(l10,1.0E-6);
                    cserrors = cserrors||ap::fp_greater(l21/l11,1.2)&&ap::fp_greater(l11,1.0E-6);
                    cserrors = cserrors||ap::fp_greater(l22/l12,1.2)&&ap::fp_greater(l12,1.0E-6);
                }
            }
            
            //
            // Build Hermite spline.
            // Test for interpolation scheme properties:
            // * values and derivatives at nodes
            // * continuous function
            // * continuous first derivative
            //
            spline1dbuildhermite(x, y, d, n, c);
            err = 0;
            for(i = 0; i <= n-1; i++)
            {
                err = ap::maxreal(err, fabs(y(i)-spline1dcalc(c, x(i))));
            }
            hserrors = hserrors||ap::fp_greater(err,threshold);
            err = 0;
            for(i = 0; i <= n-1; i++)
            {
                t = (spline1dcalc(c, x(i)+h)-spline1dcalc(c, x(i)-h))/(2*h);
                err = ap::maxreal(err, fabs(d(i)-t));
            }
            hserrors = hserrors||ap::fp_greater(err,1.0E-3);
            lconst(a, b, c, lstep, l10, l11, l12);
            lconst(a, b, c, lstep/3, l20, l21, l22);
            hserrors = hserrors||ap::fp_greater(l20/l10,1.2);
            hserrors = hserrors||ap::fp_greater(l21/l11,1.2);
            
            //
            // Build Akima spline
            // Test for general interpolation scheme properties:
            // * values at nodes
            // * continuous function
            // * continuous first derivative
            // Test for specific properties is implemented below.
            //
            if( n>=5 )
            {
                spline1dbuildakima(x, y, n, c);
                err = 0;
                for(i = 0; i <= n-1; i++)
                {
                    err = ap::maxreal(err, fabs(y(i)-spline1dcalc(c, x(i))));
                }
                aserrors = aserrors||ap::fp_greater(err,threshold);
                lconst(a, b, c, lstep, l10, l11, l12);
                lconst(a, b, c, lstep/3, l20, l21, l22);
                hserrors = hserrors||ap::fp_greater(l20/l10,1.2);
                hserrors = hserrors||ap::fp_greater(l21/l11,1.2);
            }
        }
    }
    
    //
    // Special linear spline test:
    // test for linearity between x[i] and x[i+1]
    //
    for(n = 2; n <= maxn; n++)
    {
        x.setbounds(0, n-1);
        y.setbounds(0, n-1);
        
        //
        // Prepare task
        //
        a = -1;
        b = +1;
        for(i = 0; i <= n-1; i++)
        {
            x(i) = a+(b-a)*i/(n-1);
            y(i) = 2*ap::randomreal()-1;
        }
        spline1dbuildlinear(x, y, n, c);
        
        //
        // Test
        //
        err = 0;
        for(k = 0; k <= n-2; k++)
        {
            a = x(k);
            b = x(k+1);
            for(pass = 1; pass <= passcount; pass++)
            {
                t = a+(b-a)*ap::randomreal();
                v = y(k)+(t-a)/(b-a)*(y(k+1)-y(k));
                err = ap::maxreal(err, fabs(spline1dcalc(c, t)-v));
            }
        }
        lserrors = lserrors||ap::fp_greater(err,threshold);
    }
    
    //
    // Special Akima test: test outlier sensitivity
    // Spline value at (x[i], x[i+1]) should depend from
    // f[i-2], f[i-1], f[i], f[i+1], f[i+2], f[i+3] only.
    //
    for(n = 5; n <= maxn; n++)
    {
        x.setbounds(0, n-1);
        y.setbounds(0, n-1);
        y2.setbounds(0, n-1);
        
        //
        // Prepare unperturbed Akima spline
        //
        a = -1;
        b = +1;
        for(i = 0; i <= n-1; i++)
        {
            x(i) = a+(b-a)*i/(n-1);
            y(i) = cos(1.3*ap::pi()*x(i)+0.4);
        }
        spline1dbuildakima(x, y, n, c);
        
        //
        // Process perturbed tasks
        //
        err = 0;
        for(k = 0; k <= n-1; k++)
        {
            ap::vmove(&y2(0), 1, &y(0), 1, ap::vlen(0,n-1));
            y2(k) = 5;
            spline1dbuildakima(x, y2, n, c2);
            
            //
            // Test left part independence
            //
            if( k-3>=1 )
            {
                a = -1;
                b = x(k-3);
                for(pass = 1; pass <= passcount; pass++)
                {
                    t = a+(b-a)*ap::randomreal();
                    err = ap::maxreal(err, fabs(spline1dcalc(c, t)-spline1dcalc(c2, t)));
                }
            }
            
            //
            // Test right part independence
            //
            if( k+3<=n-2 )
            {
                a = x(k+3);
                b = +1;
                for(pass = 1; pass <= passcount; pass++)
                {
                    t = a+(b-a)*ap::randomreal();
                    err = ap::maxreal(err, fabs(spline1dcalc(c, t)-spline1dcalc(c2, t)));
                }
            }
        }
        aserrors = aserrors||ap::fp_greater(err,threshold);
    }
    
    //
    // Differentiation, copy/serialize/unpack test
    //
    for(n = 2; n <= maxn; n++)
    {
        x.setbounds(0, n-1);
        y.setbounds(0, n-1);
        
        //
        // Prepare cubic spline
        //
        a = -1-ap::randomreal();
        b = +1+ap::randomreal();
        for(i = 0; i <= n-1; i++)
        {
            x(i) = a+(b-a)*i/(n-1);
            y(i) = cos(1.3*ap::pi()*x(i)+0.4);
        }
        spline1dbuildcubic(x, y, n, 2, 0.0, 2, 0.0, c);
        
        //
        // Test diff
        //
        err = 0;
        for(pass = 1; pass <= passcount; pass++)
        {
            t = a+(b-a)*ap::randomreal();
            spline1ddiff(c, t, s, ds, d2s);
            vl = spline1dcalc(c, t-h);
            vm = spline1dcalc(c, t);
            vr = spline1dcalc(c, t+h);
            err = ap::maxreal(err, fabs(s-vm));
            err = ap::maxreal(err, fabs(ds-(vr-vl)/(2*h)));
            err = ap::maxreal(err, fabs(d2s-(vr-2*vm+vl)/ap::sqr(h)));
        }
        dserrors = dserrors||ap::fp_greater(err,0.001);
        
        //
        // Test copy
        //
        unsetspline1d(c2);
        spline1dcopy(c, c2);
        err = 0;
        for(pass = 1; pass <= passcount; pass++)
        {
            t = a+(b-a)*ap::randomreal();
            err = ap::maxreal(err, fabs(spline1dcalc(c, t)-spline1dcalc(c2, t)));
        }
        cperrors = cperrors||ap::fp_greater(err,threshold);
        
        //
        // Test serialize/deserialize
        //
        unsetspline1d(c2);
        unset1d(ra2);
        spline1dserialize(c, ra, ralen);
        ra2.setlength(ralen);
        ap::vmove(&ra2(0), 1, &ra(0), 1, ap::vlen(0,ralen-1));
        spline1dunserialize(ra2, c2);
        err = 0;
        for(pass = 1; pass <= passcount; pass++)
        {
            t = a+(b-a)*ap::randomreal();
            err = ap::maxreal(err, fabs(spline1dcalc(c, t)-spline1dcalc(c2, t)));
        }
        cperrors = cperrors||ap::fp_greater(err,threshold);
        
        //
        // Test unpack
        //
        uperrors = uperrors||!testunpack(c, x);
        
        //
        // Test lin.trans.
        //
        err = 0;
        for(pass = 1; pass <= passcount; pass++)
        {
            
            //
            // LinTransX, general A
            //
            sa = 4*ap::randomreal()-2;
            sb = 2*ap::randomreal()-1;
            t = a+(b-a)*ap::randomreal();
            spline1dcopy(c, c2);
            spline1dlintransx(c2, sa, sb);
            err = ap::maxreal(err, fabs(spline1dcalc(c, t)-spline1dcalc(c2, (t-sb)/sa)));
            
            //
            // LinTransX, special case: A=0
            //
            sb = 2*ap::randomreal()-1;
            t = a+(b-a)*ap::randomreal();
            spline1dcopy(c, c2);
            spline1dlintransx(c2, double(0), sb);
            err = ap::maxreal(err, fabs(spline1dcalc(c, sb)-spline1dcalc(c2, t)));
            
            //
            // LinTransY
            //
            sa = 2*ap::randomreal()-1;
            sb = 2*ap::randomreal()-1;
            t = a+(b-a)*ap::randomreal();
            spline1dcopy(c, c2);
            spline1dlintransy(c2, sa, sb);
            err = ap::maxreal(err, fabs(sa*spline1dcalc(c, t)+sb-spline1dcalc(c2, t)));
        }
        lterrors = lterrors||ap::fp_greater(err,threshold);
    }
    
    //
    // Testing integration. Two tests are performed:
    // * approximate test (well behaved smooth function, many points,
    //   integration inside [a,b])
    // * exact test (integration of parabola, outside of [a,b]
    //
    err = 0;
    for(n = 20; n <= 35; n++)
    {
        x.setbounds(0, n-1);
        y.setbounds(0, n-1);
        for(pass = 1; pass <= passcount; pass++)
        {
            
            //
            // Prepare cubic spline
            //
            a = -1-0.2*ap::randomreal();
            b = +1+0.2*ap::randomreal();
            for(i = 0; i <= n-1; i++)
            {
                x(i) = a+(b-a)*i/(n-1);
                y(i) = sin(ap::pi()*x(i)+0.4)+exp(x(i));
            }
            bl = ap::pi()*cos(ap::pi()*a+0.4)+exp(a);
            br = ap::pi()*cos(ap::pi()*b+0.4)+exp(b);
            spline1dbuildcubic(x, y, n, 1, bl, 1, br, c);
            
            //
            // Test
            //
            t = a+(b-a)*ap::randomreal();
            v = -cos(ap::pi()*a+0.4)/ap::pi()+exp(a);
            v = -cos(ap::pi()*t+0.4)/ap::pi()+exp(t)-v;
            v = v-spline1dintegrate(c, t);
            err = ap::maxreal(err, fabs(v));
        }
    }
    ierrors = ierrors||ap::fp_greater(err,0.001);
    p0 = 2*ap::randomreal()-1;
    p1 = 2*ap::randomreal()-1;
    p2 = 2*ap::randomreal()-1;
    a = -ap::randomreal()-0.5;
    b = +ap::randomreal()+0.5;
    n = 2;
    x.setlength(n);
    y.setlength(n);
    d.setlength(n);
    x(0) = a;
    y(0) = p0+p1*a+p2*ap::sqr(a);
    d(0) = p1+2*p2*a;
    x(1) = b;
    y(1) = p0+p1*b+p2*ap::sqr(b);
    d(1) = p1+2*p2*b;
    spline1dbuildhermite(x, y, d, n, c);
    bl = ap::minreal(a, b)-fabs(b-a);
    br = ap::minreal(a, b)+fabs(b-a);
    err = 0;
    for(pass = 1; pass <= 100; pass++)
    {
        t = bl+(br-bl)*ap::randomreal();
        v = p0*t+p1*ap::sqr(t)/2+p2*ap::sqr(t)*t/3-(p0*a+p1*ap::sqr(a)/2+p2*ap::sqr(a)*a/3);
        v = v-spline1dintegrate(c, t);
        err = ap::maxreal(err, fabs(v));
    }
    ierrors = ierrors||ap::fp_greater(err,threshold);
    
    //
    // Test fitting.
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        
        //
        // Cubic splines
        // Ability to handle boundary constraints (1-4 constraints on F, dF/dx).
        //
        for(m = 4; m <= 8; m++)
        {
            for(k = 1; k <= 4; k++)
            {
                if( k>=m )
                {
                    continue;
                }
                n = 100;
                x.setlength(n);
                y.setlength(n);
                w.setlength(n);
                xc.setlength(4);
                yc.setlength(4);
                dc.setlength(4);
                sa = 1+ap::randomreal();
                sb = 2*ap::randomreal()-1;
                for(i = 0; i <= n-1; i++)
                {
                    x(i) = sa*ap::randomreal()+sb;
                    y(i) = 2*ap::randomreal()-1;
                    w(i) = 1+ap::randomreal();
                }
                xc(0) = sb;
                yc(0) = 2*ap::randomreal()-1;
                dc(0) = 0;
                xc(1) = sb;
                yc(1) = 2*ap::randomreal()-1;
                dc(1) = 1;
                xc(2) = sa+sb;
                yc(2) = 2*ap::randomreal()-1;
                dc(2) = 0;
                xc(3) = sa+sb;
                yc(3) = 2*ap::randomreal()-1;
                dc(3) = 1;
                spline1dfitcubicwc(x, y, w, n, xc, yc, dc, k, m, info, c, rep);
                if( info<=0 )
                {
                    fiterrors = true;
                }
                else
                {
                    
                    //
                    // Check that constraints are satisfied
                    //
                    for(i = 0; i <= k-1; i++)
                    {
                        spline1ddiff(c, xc(i), s, ds, d2s);
                        if( dc(i)==0 )
                        {
                            fiterrors = fiterrors||ap::fp_greater(fabs(s-yc(i)),threshold);
                        }
                        if( dc(i)==1 )
                        {
                            fiterrors = fiterrors||ap::fp_greater(fabs(ds-yc(i)),threshold);
                        }
                        if( dc(i)==2 )
                        {
                            fiterrors = fiterrors||ap::fp_greater(fabs(d2s-yc(i)),threshold);
                        }
                    }
                }
            }
        }
        
        //
        // Cubic splines
        // Ability to handle one internal constraint
        //
        for(m = 4; m <= 8; m++)
        {
            n = 100;
            x.setlength(n);
            y.setlength(n);
            w.setlength(n);
            xc.setlength(1);
            yc.setlength(1);
            dc.setlength(1);
            sa = 1+ap::randomreal();
            sb = 2*ap::randomreal()-1;
            for(i = 0; i <= n-1; i++)
            {
                x(i) = sa*ap::randomreal()+sb;
                y(i) = 2*ap::randomreal()-1;
                w(i) = 1+ap::randomreal();
            }
            xc(0) = sa*ap::randomreal()+sb;
            yc(0) = 2*ap::randomreal()-1;
            dc(0) = ap::randominteger(2);
            spline1dfitcubicwc(x, y, w, n, xc, yc, dc, 1, m, info, c, rep);
            if( info<=0 )
            {
                fiterrors = true;
            }
            else
            {
                
                //
                // Check that constraints are satisfied
                //
                spline1ddiff(c, xc(0), s, ds, d2s);
                if( dc(0)==0 )
                {
                    fiterrors = fiterrors||ap::fp_greater(fabs(s-yc(0)),threshold);
                }
                if( dc(0)==1 )
                {
                    fiterrors = fiterrors||ap::fp_greater(fabs(ds-yc(0)),threshold);
                }
                if( dc(0)==2 )
                {
                    fiterrors = fiterrors||ap::fp_greater(fabs(d2s-yc(0)),threshold);
                }
            }
        }
        
        //
        // Hermite splines
        // Ability to handle boundary constraints (1-4 constraints on F, dF/dx).
        //
        for(m = 4; m <= 8; m++)
        {
            for(k = 1; k <= 4; k++)
            {
                if( k>=m )
                {
                    continue;
                }
                if( m%2!=0 )
                {
                    continue;
                }
                n = 100;
                x.setlength(n);
                y.setlength(n);
                w.setlength(n);
                xc.setlength(4);
                yc.setlength(4);
                dc.setlength(4);
                sa = 1+ap::randomreal();
                sb = 2*ap::randomreal()-1;
                for(i = 0; i <= n-1; i++)
                {
                    x(i) = sa*ap::randomreal()+sb;
                    y(i) = 2*ap::randomreal()-1;
                    w(i) = 1+ap::randomreal();
                }
                xc(0) = sb;
                yc(0) = 2*ap::randomreal()-1;
                dc(0) = 0;
                xc(1) = sb;
                yc(1) = 2*ap::randomreal()-1;
                dc(1) = 1;
                xc(2) = sa+sb;
                yc(2) = 2*ap::randomreal()-1;
                dc(2) = 0;
                xc(3) = sa+sb;
                yc(3) = 2*ap::randomreal()-1;
                dc(3) = 1;
                spline1dfithermitewc(x, y, w, n, xc, yc, dc, k, m, info, c, rep);
                if( info<=0 )
                {
                    fiterrors = true;
                }
                else
                {
                    
                    //
                    // Check that constraints are satisfied
                    //
                    for(i = 0; i <= k-1; i++)
                    {
                        spline1ddiff(c, xc(i), s, ds, d2s);
                        if( dc(i)==0 )
                        {
                            fiterrors = fiterrors||ap::fp_greater(fabs(s-yc(i)),threshold);
                        }
                        if( dc(i)==1 )
                        {
                            fiterrors = fiterrors||ap::fp_greater(fabs(ds-yc(i)),threshold);
                        }
                        if( dc(i)==2 )
                        {
                            fiterrors = fiterrors||ap::fp_greater(fabs(d2s-yc(i)),threshold);
                        }
                    }
                }
            }
        }
        
        //
        // Hermite splines
        // Ability to handle one internal constraint
        //
        for(m = 4; m <= 8; m++)
        {
            if( m%2!=0 )
            {
                continue;
            }
            n = 100;
            x.setlength(n);
            y.setlength(n);
            w.setlength(n);
            xc.setlength(1);
            yc.setlength(1);
            dc.setlength(1);
            sa = 1+ap::randomreal();
            sb = 2*ap::randomreal()-1;
            for(i = 0; i <= n-1; i++)
            {
                x(i) = sa*ap::randomreal()+sb;
                y(i) = 2*ap::randomreal()-1;
                w(i) = 1+ap::randomreal();
            }
            xc(0) = sa*ap::randomreal()+sb;
            yc(0) = 2*ap::randomreal()-1;
            dc(0) = ap::randominteger(2);
            spline1dfithermitewc(x, y, w, n, xc, yc, dc, 1, m, info, c, rep);
            if( info<=0 )
            {
                fiterrors = true;
            }
            else
            {
                
                //
                // Check that constraints are satisfied
                //
                spline1ddiff(c, xc(0), s, ds, d2s);
                if( dc(0)==0 )
                {
                    fiterrors = fiterrors||ap::fp_greater(fabs(s-yc(0)),threshold);
                }
                if( dc(0)==1 )
                {
                    fiterrors = fiterrors||ap::fp_greater(fabs(ds-yc(0)),threshold);
                }
                if( dc(0)==2 )
                {
                    fiterrors = fiterrors||ap::fp_greater(fabs(d2s-yc(0)),threshold);
                }
            }
        }
    }
    for(m = 4; m <= 8; m++)
    {
        for(stype = 0; stype <= 1; stype++)
        {
            for(pass = 1; pass <= passcount; pass++)
            {
                if( stype==1&&m%2!=0 )
                {
                    continue;
                }
                
                //
                // cubic/Hermite spline fitting:
                // * generate "template spline" C2
                // * generate 2*N points from C2, such that result of
                //   ideal fit should be equal to C2
                // * fit, store in C
                // * compare C and C2
                //
                sa = 1+ap::randomreal();
                sb = 2*ap::randomreal()-1;
                if( stype==0 )
                {
                    x.setlength(m-2);
                    y.setlength(m-2);
                    for(i = 0; i <= m-2-1; i++)
                    {
                        x(i) = sa*i/(m-2-1)+sb;
                        y(i) = 2*ap::randomreal()-1;
                    }
                    spline1dbuildcubic(x, y, m-2, 1, 2*ap::randomreal()-1, 1, 2*ap::randomreal()-1, c2);
                }
                if( stype==1 )
                {
                    x.setlength(m/2);
                    y.setlength(m/2);
                    d.setlength(m/2);
                    for(i = 0; i <= m/2-1; i++)
                    {
                        x(i) = sa*i/(m/2-1)+sb;
                        y(i) = 2*ap::randomreal()-1;
                        d(i) = 2*ap::randomreal()-1;
                    }
                    spline1dbuildhermite(x, y, d, m/2, c2);
                }
                n = 50;
                x.setlength(2*n);
                y.setlength(2*n);
                w.setlength(2*n);
                for(i = 0; i <= n-1; i++)
                {
                    
                    //
                    // "if i=0" and "if i=1" are needed to
                    // synchronize interval size for C2 and
                    // spline being fitted (i.e. C).
                    //
                    t = ap::randomreal();
                    x(i) = sa*ap::randomreal()+sb;
                    if( i==0 )
                    {
                        x(i) = sb;
                    }
                    if( i==1 )
                    {
                        x(i) = sa+sb;
                    }
                    v = spline1dcalc(c2, x(i));
                    y(i) = v+t;
                    w(i) = 1+ap::randomreal();
                    x(n+i) = x(i);
                    y(n+i) = v-t;
                    w(n+i) = w(i);
                }
                if( stype==0 )
                {
                    spline1dfitcubicwc(x, y, w, 2*n, xc, yc, dc, 0, m, info, c, rep);
                }
                if( stype==1 )
                {
                    spline1dfithermitewc(x, y, w, 2*n, xc, yc, dc, 0, m, info, c, rep);
                }
                if( info<=0 )
                {
                    fiterrors = true;
                }
                else
                {
                    for(i = 0; i <= n-1; i++)
                    {
                        v = sa*ap::randomreal()+sb;
                        fiterrors = fiterrors||ap::fp_greater(fabs(spline1dcalc(c, v)-spline1dcalc(c2, v)),threshold);
                    }
                }
            }
        }
    }
    for(m = 4; m <= 8; m++)
    {
        for(pass = 1; pass <= passcount; pass++)
        {
            
            //
            // prepare points/weights
            //
            sa = 1+ap::randomreal();
            sb = 2*ap::randomreal()-1;
            n = 10+ap::randominteger(10);
            x.setlength(n);
            y.setlength(n);
            w.setlength(n);
            for(i = 0; i <= n-1; i++)
            {
                x(i) = sa*ap::randomreal()+sb;
                y(i) = 2*ap::randomreal()-1;
                w(i) = 1;
            }
            
            //
            // Fit cubic with unity weights, without weights, then compare
            //
            if( m>=4 )
            {
                spline1dfitcubicwc(x, y, w, n, xc, yc, dc, 0, m, info1, c, rep);
                spline1dfitcubic(x, y, n, m, info2, c2, rep2);
                if( info1<=0||info2<=0 )
                {
                    fiterrors = true;
                }
                else
                {
                    for(i = 0; i <= n-1; i++)
                    {
                        v = sa*ap::randomreal()+sb;
                        fiterrors = fiterrors||ap::fp_neq(spline1dcalc(c, v),spline1dcalc(c2, v));
                        fiterrors = fiterrors||ap::fp_neq(rep.taskrcond,rep2.taskrcond);
                        fiterrors = fiterrors||ap::fp_neq(rep.rmserror,rep2.rmserror);
                        fiterrors = fiterrors||ap::fp_neq(rep.avgerror,rep2.avgerror);
                        fiterrors = fiterrors||ap::fp_neq(rep.avgrelerror,rep2.avgrelerror);
                        fiterrors = fiterrors||ap::fp_neq(rep.maxerror,rep2.maxerror);
                    }
                }
            }
            
            //
            // Fit Hermite with unity weights, without weights, then compare
            //
            if( m>=4&&m%2==0 )
            {
                spline1dfithermitewc(x, y, w, n, xc, yc, dc, 0, m, info1, c, rep);
                spline1dfithermite(x, y, n, m, info2, c2, rep2);
                if( info1<=0||info2<=0 )
                {
                    fiterrors = true;
                }
                else
                {
                    for(i = 0; i <= n-1; i++)
                    {
                        v = sa*ap::randomreal()+sb;
                        fiterrors = fiterrors||ap::fp_neq(spline1dcalc(c, v),spline1dcalc(c2, v));
                        fiterrors = fiterrors||ap::fp_neq(rep.taskrcond,rep2.taskrcond);
                        fiterrors = fiterrors||ap::fp_neq(rep.rmserror,rep2.rmserror);
                        fiterrors = fiterrors||ap::fp_neq(rep.avgerror,rep2.avgerror);
                        fiterrors = fiterrors||ap::fp_neq(rep.avgrelerror,rep2.avgrelerror);
                        fiterrors = fiterrors||ap::fp_neq(rep.maxerror,rep2.maxerror);
                    }
                }
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
        // Test cubic fitting
        //
        spline1dfitcubic(x, y, 4, 4, info, c, rep);
        if( info<=0 )
        {
            fiterrors = true;
        }
        else
        {
            s = spline1dcalc(c, double(0));
            fiterrors = fiterrors||ap::fp_greater(fabs(s-v),threshold);
            fiterrors = fiterrors||ap::fp_greater(fabs(rep.rmserror-refrms),threshold);
            fiterrors = fiterrors||ap::fp_greater(fabs(rep.avgerror-refavg),threshold);
            fiterrors = fiterrors||ap::fp_greater(fabs(rep.avgrelerror-refavgrel),threshold);
            fiterrors = fiterrors||ap::fp_greater(fabs(rep.maxerror-refmax),threshold);
        }
        
        //
        // Test cubic fitting
        //
        spline1dfithermite(x, y, 4, 4, info, c, rep);
        if( info<=0 )
        {
            fiterrors = true;
        }
        else
        {
            s = spline1dcalc(c, double(0));
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
    waserrors = lserrors||cserrors||hserrors||aserrors||dserrors||cperrors||uperrors||lterrors||ierrors||fiterrors;
    if( !silent )
    {
        printf("TESTING SPLINE INTERPOLATION\n");
        
        //
        // Normal tests
        //
        printf("LINEAR SPLINE TEST:                      ");
        if( lserrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("CUBIC SPLINE TEST:                       ");
        if( cserrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("HERMITE SPLINE TEST:                     ");
        if( hserrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("AKIMA SPLINE TEST:                       ");
        if( aserrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("DIFFERENTIATION TEST:                    ");
        if( dserrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("COPY/SERIALIZATION TEST:                 ");
        if( cperrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("UNPACK TEST:                             ");
        if( uperrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("LIN.TRANS. TEST:                         ");
        if( lterrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("INTEGRATION TEST:                        ");
        if( ierrors )
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


/*************************************************************************
Lipschitz constants for spline inself, first and second derivatives.
*************************************************************************/
static void lconst(double a,
     double b,
     const spline1dinterpolant& c,
     double lstep,
     double& l0,
     double& l1,
     double& l2)
{
    double t;
    double vl;
    double vm;
    double vr;
    double prevf;
    double prevd;
    double prevd2;
    double f;
    double d;
    double d2;

    l0 = 0;
    l1 = 0;
    l2 = 0;
    t = a-0.1;
    vl = spline1dcalc(c, t-2*lstep);
    vm = spline1dcalc(c, t-lstep);
    vr = spline1dcalc(c, t);
    f = vm;
    d = (vr-vl)/(2*lstep);
    d2 = (vr-2*vm+vl)/ap::sqr(lstep);
    while(ap::fp_less_eq(t,b+0.1))
    {
        prevf = f;
        prevd = d;
        prevd2 = d2;
        vl = vm;
        vm = vr;
        vr = spline1dcalc(c, t+lstep);
        f = vm;
        d = (vr-vl)/(2*lstep);
        d2 = (vr-2*vm+vl)/ap::sqr(lstep);
        l0 = ap::maxreal(l0, fabs((f-prevf)/lstep));
        l1 = ap::maxreal(l1, fabs((d-prevd)/lstep));
        l2 = ap::maxreal(l2, fabs((d2-prevd2)/lstep));
        t = t+lstep;
    }
}


/*************************************************************************
Unpack testing
*************************************************************************/
static bool testunpack(const spline1dinterpolant& c,
     const ap::real_1d_array& x)
{
    bool result;
    int i;
    int n;
    double err;
    double t;
    double v1;
    double v2;
    int pass;
    int passcount;
    ap::real_2d_array tbl;

    passcount = 20;
    err = 0;
    spline1dunpack(c, n, tbl);
    for(i = 0; i <= n-2; i++)
    {
        for(pass = 1; pass <= passcount; pass++)
        {
            t = ap::randomreal()*(tbl(i,1)-tbl(i,0));
            v1 = tbl(i,2)+t*tbl(i,3)+ap::sqr(t)*tbl(i,4)+t*ap::sqr(t)*tbl(i,5);
            v2 = spline1dcalc(c, tbl(i,0)+t);
            err = ap::maxreal(err, fabs(v1-v2));
        }
    }
    for(i = 0; i <= n-2; i++)
    {
        err = ap::maxreal(err, fabs(x(i)-tbl(i,0)));
    }
    for(i = 0; i <= n-2; i++)
    {
        err = ap::maxreal(err, fabs(x(i+1)-tbl(i,1)));
    }
    result = ap::fp_less(err,100*ap::machineepsilon);
    return result;
}


/*************************************************************************
Unset spline, i.e. initialize it with random garbage
*************************************************************************/
static void unsetspline1d(spline1dinterpolant& c)
{
    ap::real_1d_array x;
    ap::real_1d_array y;
    ap::real_1d_array d;

    x.setlength(2);
    y.setlength(2);
    d.setlength(2);
    x(0) = -1;
    y(0) = ap::randomreal();
    d(0) = ap::randomreal();
    x(1) = 1;
    y(1) = ap::randomreal();
    d(1) = ap::randomreal();
    spline1dbuildhermite(x, y, d, 2, c);
}


/*************************************************************************
Unsets real vector
*************************************************************************/
static void unset1d(ap::real_1d_array& x)
{

    x.setlength(1);
    x(0) = 2*ap::randomreal()-1;
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
bool testspline1dunit_test_silent()
{
    bool result;

    result = testsplineinterpolation(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testspline1dunit_test()
{
    bool result;

    result = testsplineinterpolation(false);
    return result;
}






 
#include <stdio.h>
#include "llstestunit.h"

static bool isglssolution(int n,
     int m,
     int k,
     const ap::real_1d_array& y,
     const ap::real_1d_array& w,
     const ap::real_2d_array& fmatrix,
     const ap::real_2d_array& cmatrix,
     ap::real_1d_array c);
static double getglserror(int n,
     int m,
     const ap::real_1d_array& y,
     const ap::real_1d_array& w,
     const ap::real_2d_array& fmatrix,
     const ap::real_1d_array& c);
static void fitlinearnonlinear(int m,
     bool gradonly,
     const ap::real_2d_array& xy,
     lsfitstate& state,
     bool& nlserrors);

bool testlls(bool silent)
{
    bool result;
    bool waserrors;
    bool llserrors;
    bool nlserrors;
    double threshold;
    double nlthreshold;
    int maxn;
    int maxm;
    int passcount;
    int n;
    int m;
    int i;
    int j;
    int k;
    int pass;
    double xscale;
    ap::real_1d_array x;
    ap::real_1d_array y;
    ap::real_1d_array w;
    ap::real_1d_array w2;
    ap::real_1d_array c;
    ap::real_1d_array c2;
    ap::real_2d_array a;
    ap::real_2d_array a2;
    ap::real_2d_array cm;
    double v;
    double v1;
    double v2;
    lsfitreport rep;
    lsfitreport rep2;
    int info;
    int info2;
    double refrms;
    double refavg;
    double refavgrel;
    double refmax;
    lsfitstate state;

    waserrors = false;
    llserrors = false;
    nlserrors = false;
    threshold = 10000*ap::machineepsilon;
    nlthreshold = 0.00001;
    maxn = 6;
    maxm = 6;
    passcount = 4;
    
    //
    // Testing unconstrained least squares (linear/nonlinear)
    //
    for(n = 1; n <= maxn; n++)
    {
        for(m = 1; m <= maxm; m++)
        {
            for(pass = 1; pass <= passcount; pass++)
            {
                
                //
                // Solve non-degenerate linear least squares task
                // Use Chebyshev basis. Its condition number is very good.
                //
                a.setlength(n, m);
                x.setlength(n);
                y.setlength(n);
                w.setlength(n);
                xscale = 0.9+0.1*ap::randomreal();
                for(i = 0; i <= n-1; i++)
                {
                    if( n==1 )
                    {
                        x(i) = 2*ap::randomreal()-1;
                    }
                    else
                    {
                        x(i) = xscale*(double(2*i)/double(n-1)-1);
                    }
                    y(i) = 3*x(i)+exp(x(i));
                    w(i) = 1+ap::randomreal();
                    a(i,0) = 1;
                    if( m>1 )
                    {
                        a(i,1) = x(i);
                    }
                    for(j = 2; j <= m-1; j++)
                    {
                        a(i,j) = 2*x(i)*a(i,j-1)-a(i,j-2);
                    }
                }
                
                //
                // 1. test weighted fitting (optimality)
                // 2. Solve degenerate least squares task built on the basis
                //    of previous task
                //
                lsfitlinearw(y, w, a, n, m, info, c, rep);
                if( info<=0 )
                {
                    llserrors = true;
                }
                else
                {
                    llserrors = llserrors||!isglssolution(n, m, 0, y, w, a, cm, c);
                }
                a2.setlength(n, 2*m);
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= m-1; j++)
                    {
                        a2(i,2*j+0) = a(i,j);
                        a2(i,2*j+1) = a(i,j);
                    }
                }
                lsfitlinearw(y, w, a2, n, 2*m, info, c2, rep);
                if( info<=0 )
                {
                    llserrors = true;
                }
                else
                {
                    
                    //
                    // test answer correctness using design matrix properties
                    // and previous task solution
                    //
                    for(j = 0; j <= m-1; j++)
                    {
                        llserrors = llserrors||ap::fp_greater(fabs(c2(2*j+0)+c2(2*j+1)-c(j)),threshold);
                    }
                }
                
                //
                // test non-weighted fitting
                //
                w2.setlength(n);
                for(i = 0; i <= n-1; i++)
                {
                    w2(i) = 1;
                }
                lsfitlinearw(y, w2, a, n, m, info, c, rep);
                lsfitlinear(y, a, n, m, info2, c2, rep2);
                if( info<=0||info2<=0 )
                {
                    llserrors = true;
                }
                else
                {
                    
                    //
                    // test answer correctness
                    //
                    for(j = 0; j <= m-1; j++)
                    {
                        llserrors = llserrors||ap::fp_greater(fabs(c(j)-c2(j)),threshold);
                    }
                    llserrors = llserrors||ap::fp_greater(fabs(rep.taskrcond-rep2.taskrcond),threshold);
                }
                
                //
                // test nonlinear fitting on the linear task
                // (only non-degenerate task are tested)
                // and compare with answer from linear fitting subroutine
                //
                if( n>=m )
                {
                    c2.setlength(m);
                    
                    //
                    // test gradient-only or Hessian-based weighted fitting
                    //
                    lsfitlinearw(y, w, a, n, m, info, c, rep);
                    for(i = 0; i <= m-1; i++)
                    {
                        c2(i) = 2*ap::randomreal()-1;
                    }
                    lsfitnonlinearwfg(a, y, w, c2, n, m, m, ap::fp_greater(ap::randomreal(),0.5), state);
                    lsfitnonlinearsetcond(state, 0.0, nlthreshold, 0);
                    fitlinearnonlinear(m, true, a, state, nlserrors);
                    lsfitnonlinearresults(state, info, c2, rep2);
                    if( info<=0 )
                    {
                        nlserrors = true;
                    }
                    else
                    {
                        for(i = 0; i <= m-1; i++)
                        {
                            nlserrors = nlserrors||ap::fp_greater(fabs(c(i)-c2(i)),100*nlthreshold);
                        }
                    }
                    for(i = 0; i <= m-1; i++)
                    {
                        c2(i) = 2*ap::randomreal()-1;
                    }
                    lsfitnonlinearwfgh(a, y, w, c2, n, m, m, state);
                    lsfitnonlinearsetcond(state, 0.0, nlthreshold, 0);
                    fitlinearnonlinear(m, false, a, state, nlserrors);
                    lsfitnonlinearresults(state, info, c2, rep2);
                    if( info<=0 )
                    {
                        nlserrors = true;
                    }
                    else
                    {
                        for(i = 0; i <= m-1; i++)
                        {
                            nlserrors = nlserrors||ap::fp_greater(fabs(c(i)-c2(i)),100*nlthreshold);
                        }
                    }
                    
                    //
                    // test gradient-only or Hessian-based fitting without weights
                    //
                    lsfitlinear(y, a, n, m, info, c, rep);
                    for(i = 0; i <= m-1; i++)
                    {
                        c2(i) = 2*ap::randomreal()-1;
                    }
                    lsfitnonlinearfg(a, y, c2, n, m, m, ap::fp_greater(ap::randomreal(),0.5), state);
                    lsfitnonlinearsetcond(state, 0.0, nlthreshold, 0);
                    fitlinearnonlinear(m, true, a, state, nlserrors);
                    lsfitnonlinearresults(state, info, c2, rep2);
                    if( info<=0 )
                    {
                        nlserrors = true;
                    }
                    else
                    {
                        for(i = 0; i <= m-1; i++)
                        {
                            nlserrors = nlserrors||ap::fp_greater(fabs(c(i)-c2(i)),100*nlthreshold);
                        }
                    }
                    for(i = 0; i <= m-1; i++)
                    {
                        c2(i) = 2*ap::randomreal()-1;
                    }
                    lsfitnonlinearfgh(a, y, c2, n, m, m, state);
                    lsfitnonlinearsetcond(state, 0.0, nlthreshold, 0);
                    fitlinearnonlinear(m, false, a, state, nlserrors);
                    lsfitnonlinearresults(state, info, c2, rep2);
                    if( info<=0 )
                    {
                        nlserrors = true;
                    }
                    else
                    {
                        for(i = 0; i <= m-1; i++)
                        {
                            nlserrors = nlserrors||ap::fp_greater(fabs(c(i)-c2(i)),100*nlthreshold);
                        }
                    }
                }
            }
        }
        
        //
        // test correctness of the RCond field
        //
        a.setbounds(0, n-1, 0, n-1);
        x.setbounds(0, n-1);
        y.setbounds(0, n-1);
        w.setbounds(0, n-1);
        v1 = ap::maxrealnumber;
        v2 = ap::minrealnumber;
        for(i = 0; i <= n-1; i++)
        {
            x(i) = 0.1+0.9*ap::randomreal();
            y(i) = 0.1+0.9*ap::randomreal();
            w(i) = 1;
            for(j = 0; j <= n-1; j++)
            {
                if( i==j )
                {
                    a(i,i) = 0.1+0.9*ap::randomreal();
                    v1 = ap::minreal(v1, a(i,i));
                    v2 = ap::maxreal(v2, a(i,i));
                }
                else
                {
                    a(i,j) = 0;
                }
            }
        }
        lsfitlinearw(y, w, a, n, n, info, c, rep);
        if( info<=0 )
        {
            llserrors = true;
        }
        else
        {
            llserrors = llserrors||ap::fp_greater(fabs(rep.taskrcond-v1/v2),threshold);
        }
    }
    
    //
    // Test constrained least squares
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        for(n = 1; n <= maxn; n++)
        {
            for(m = 1; m <= maxm; m++)
            {
                
                //
                // test for K<>0
                //
                for(k = 1; k <= m-1; k++)
                {
                    
                    //
                    // Prepare Chebyshev basis. Its condition number is very good.
                    // Prepare constraints (random numbers)
                    //
                    a.setlength(n, m);
                    x.setlength(n);
                    y.setlength(n);
                    w.setlength(n);
                    xscale = 0.9+0.1*ap::randomreal();
                    for(i = 0; i <= n-1; i++)
                    {
                        if( n==1 )
                        {
                            x(i) = 2*ap::randomreal()-1;
                        }
                        else
                        {
                            x(i) = xscale*(double(2*i)/double(n-1)-1);
                        }
                        y(i) = 3*x(i)+exp(x(i));
                        w(i) = 1+ap::randomreal();
                        a(i,0) = 1;
                        if( m>1 )
                        {
                            a(i,1) = x(i);
                        }
                        for(j = 2; j <= m-1; j++)
                        {
                            a(i,j) = 2*x(i)*a(i,j-1)-a(i,j-2);
                        }
                    }
                    cm.setlength(k, m+1);
                    for(i = 0; i <= k-1; i++)
                    {
                        for(j = 0; j <= m; j++)
                        {
                            cm(i,j) = 2*ap::randomreal()-1;
                        }
                    }
                    
                    //
                    // Solve constrained task
                    //
                    lsfitlinearwc(y, w, a, cm, n, m, k, info, c, rep);
                    if( info<=0 )
                    {
                        llserrors = true;
                    }
                    else
                    {
                        llserrors = llserrors||!isglssolution(n, m, k, y, w, a, cm, c);
                    }
                    
                    //
                    // test non-weighted fitting
                    //
                    w2.setlength(n);
                    for(i = 0; i <= n-1; i++)
                    {
                        w2(i) = 1;
                    }
                    lsfitlinearwc(y, w2, a, cm, n, m, k, info, c, rep);
                    lsfitlinearc(y, a, cm, n, m, k, info2, c2, rep2);
                    if( info<=0||info2<=0 )
                    {
                        llserrors = true;
                    }
                    else
                    {
                        
                        //
                        // test answer correctness
                        //
                        for(j = 0; j <= m-1; j++)
                        {
                            llserrors = llserrors||ap::fp_greater(fabs(c(j)-c2(j)),threshold);
                        }
                        llserrors = llserrors||ap::fp_greater(fabs(rep.taskrcond-rep2.taskrcond),threshold);
                    }
                }
            }
        }
    }
    
    //
    // nonlinear task for nonlinear fitting:
    //
    //     f(X,C) = 1/(1+C*X^2),
    //     C(true) = 2.
    //
    n = 100;
    c.setlength(1);
    c(0) = 1+2*ap::randomreal();
    a.setlength(n, 1);
    y.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        a(i,0) = 4*ap::randomreal()-2;
        y(i) = 1/(1+2*ap::sqr(a(i,0)));
    }
    lsfitnonlinearfg(a, y, c, n, 1, 1, true, state);
    lsfitnonlinearsetcond(state, 0.0, nlthreshold, 0);
    while(lsfitnonlineariteration(state))
    {
        if( state.needf )
        {
            state.f = 1/(1+state.c(0)*ap::sqr(state.x(0)));
        }
        if( state.needfg )
        {
            state.f = 1/(1+state.c(0)*ap::sqr(state.x(0)));
            state.g(0) = -ap::sqr(state.x(0))/ap::sqr(1+state.c(0)*ap::sqr(state.x(0)));
        }
    }
    lsfitnonlinearresults(state, info, c, rep);
    if( info<=0 )
    {
        nlserrors = true;
    }
    else
    {
        nlserrors = nlserrors||ap::fp_greater(fabs(c(0)-2),100*nlthreshold);
    }
    
    //
    // solve simple task (fitting by constant function) and check
    // correctness of the errors calculated by subroutines
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        
        //
        // test on task with non-zero Yi
        //
        n = 4;
        v1 = ap::randomreal();
        v2 = ap::randomreal();
        v = 1+ap::randomreal();
        c.setlength(1);
        c(0) = 1+2*ap::randomreal();
        a.setlength(4, 1);
        y.setlength(4);
        a(0,0) = 1;
        y(0) = v-v2;
        a(1,0) = 1;
        y(1) = v-v1;
        a(2,0) = 1;
        y(2) = v+v1;
        a(3,0) = 1;
        y(3) = v+v2;
        refrms = sqrt((ap::sqr(v1)+ap::sqr(v2))/2);
        refavg = (fabs(v1)+fabs(v2))/2;
        refavgrel = 0.25*(fabs(v2)/fabs(v-v2)+fabs(v1)/fabs(v-v1)+fabs(v1)/fabs(v+v1)+fabs(v2)/fabs(v+v2));
        refmax = ap::maxreal(v1, v2);
        
        //
        // Test LLS
        //
        lsfitlinear(y, a, 4, 1, info, c, rep);
        if( info<=0 )
        {
            llserrors = true;
        }
        else
        {
            llserrors = llserrors||ap::fp_greater(fabs(c(0)-v),threshold);
            llserrors = llserrors||ap::fp_greater(fabs(rep.rmserror-refrms),threshold);
            llserrors = llserrors||ap::fp_greater(fabs(rep.avgerror-refavg),threshold);
            llserrors = llserrors||ap::fp_greater(fabs(rep.avgrelerror-refavgrel),threshold);
            llserrors = llserrors||ap::fp_greater(fabs(rep.maxerror-refmax),threshold);
        }
        
        //
        // Test NLS
        //
        lsfitnonlinearfg(a, y, c, 4, 1, 1, true, state);
        lsfitnonlinearsetcond(state, 0.0, nlthreshold, 0);
        while(lsfitnonlineariteration(state))
        {
            if( state.needf )
            {
                state.f = state.c(0);
            }
            if( state.needfg )
            {
                state.f = state.c(0);
                state.g(0) = 1;
            }
        }
        lsfitnonlinearresults(state, info, c, rep);
        if( info<=0 )
        {
            nlserrors = true;
        }
        else
        {
            nlserrors = nlserrors||ap::fp_greater(fabs(c(0)-v),threshold);
            nlserrors = nlserrors||ap::fp_greater(fabs(rep.rmserror-refrms),threshold);
            nlserrors = nlserrors||ap::fp_greater(fabs(rep.avgerror-refavg),threshold);
            nlserrors = nlserrors||ap::fp_greater(fabs(rep.avgrelerror-refavgrel),threshold);
            nlserrors = nlserrors||ap::fp_greater(fabs(rep.maxerror-refmax),threshold);
        }
    }
    
    //
    // report
    //
    waserrors = llserrors||nlserrors;
    if( !silent )
    {
        printf("TESTING LEAST SQUARES\n");
        printf("LINEAR LEAST SQUARES:                    ");
        if( llserrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("NON-LINEAR LEAST SQUARES:                ");
        if( nlserrors )
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
Tests whether C is solution of (possibly) constrained LLS problem
*************************************************************************/
static bool isglssolution(int n,
     int m,
     int k,
     const ap::real_1d_array& y,
     const ap::real_1d_array& w,
     const ap::real_2d_array& fmatrix,
     const ap::real_2d_array& cmatrix,
     ap::real_1d_array c)
{
    bool result;
    int i;
    int j;
    ap::real_1d_array c2;
    ap::real_1d_array sv;
    ap::real_1d_array deltac;
    ap::real_1d_array deltaproj;
    ap::real_2d_array u;
    ap::real_2d_array vt;
    double v;
    double s1;
    double s2;
    double s3;
    double delta;
    double threshold;

    
    //
    // Setup.
    // Threshold is small because CMatrix may be ill-conditioned
    //
    delta = 0.001;
    threshold = sqrt(ap::machineepsilon);
    c2.setlength(m);
    deltac.setlength(m);
    deltaproj.setlength(m);
    
    //
    // test whether C is feasible point or not (projC must be close to C)
    //
    for(i = 0; i <= k-1; i++)
    {
        v = ap::vdotproduct(&cmatrix(i, 0), 1, &c(0), 1, ap::vlen(0,m-1));
        if( ap::fp_greater(fabs(v-cmatrix(i,m)),threshold) )
        {
            result = false;
            return result;
        }
    }
    
    //
    // find orthogonal basis of Null(CMatrix) (stored in rows from K to M-1)
    //
    if( k>0 )
    {
        rmatrixsvd(cmatrix, k, m, 0, 2, 2, sv, u, vt);
    }
    
    //
    // Test result
    //
    result = true;
    s1 = getglserror(n, m, y, w, fmatrix, c);
    for(j = 0; j <= m-1; j++)
    {
        
        //
        // prepare modification of C which leave us in the feasible set.
        //
        // let deltaC be increment on Jth coordinate, then project
        // deltaC in the Null(CMatrix) and store result in DeltaProj
        //
        ap::vmove(&c2(0), 1, &c(0), 1, ap::vlen(0,m-1));
        for(i = 0; i <= m-1; i++)
        {
            if( i==j )
            {
                deltac(i) = delta;
            }
            else
            {
                deltac(i) = 0;
            }
        }
        if( k==0 )
        {
            ap::vmove(&deltaproj(0), 1, &deltac(0), 1, ap::vlen(0,m-1));
        }
        else
        {
            for(i = 0; i <= m-1; i++)
            {
                deltaproj(i) = 0;
            }
            for(i = k; i <= m-1; i++)
            {
                v = ap::vdotproduct(&vt(i, 0), 1, &deltac(0), 1, ap::vlen(0,m-1));
                ap::vadd(&deltaproj(0), 1, &vt(i, 0), 1, ap::vlen(0,m-1), v);
            }
        }
        
        //
        // now we have DeltaProj such that if C is feasible,
        // then C+DeltaProj is feasible too
        //
        ap::vmove(&c2(0), 1, &c(0), 1, ap::vlen(0,m-1));
        ap::vadd(&c2(0), 1, &deltaproj(0), 1, ap::vlen(0,m-1));
        s2 = getglserror(n, m, y, w, fmatrix, c2);
        ap::vmove(&c2(0), 1, &c(0), 1, ap::vlen(0,m-1));
        ap::vsub(&c2(0), 1, &deltaproj(0), 1, ap::vlen(0,m-1));
        s3 = getglserror(n, m, y, w, fmatrix, c2);
        result = result&&ap::fp_greater_eq(s2,s1/(1+threshold))&&ap::fp_greater_eq(s3,s1/(1+threshold));
    }
    return result;
}


/*************************************************************************
Tests whether C is solution of LLS problem
*************************************************************************/
static double getglserror(int n,
     int m,
     const ap::real_1d_array& y,
     const ap::real_1d_array& w,
     const ap::real_2d_array& fmatrix,
     const ap::real_1d_array& c)
{
    double result;
    int i;
    double v;

    result = 0;
    for(i = 0; i <= n-1; i++)
    {
        v = ap::vdotproduct(&fmatrix(i, 0), 1, &c(0), 1, ap::vlen(0,m-1));
        result = result+ap::sqr(w(i)*(v-y(i)));
    }
    return result;
}


/*************************************************************************
Subroutine for nonlinear fitting of linear problem
*************************************************************************/
static void fitlinearnonlinear(int m,
     bool gradonly,
     const ap::real_2d_array& xy,
     lsfitstate& state,
     bool& nlserrors)
{
    int i;
    int j;
    double v;

    while(lsfitnonlineariteration(state))
    {
        
        //
        // assume that one and only one of flags is set
        // test that we didn't request hessian in hessian-free setting
        //
        if( gradonly&&state.needfgh )
        {
            nlserrors = true;
        }
        i = 0;
        if( state.needf )
        {
            i = i+1;
        }
        if( state.needfg )
        {
            i = i+1;
        }
        if( state.needfgh )
        {
            i = i+1;
        }
        if( i!=1 )
        {
            nlserrors = true;
        }
        
        //
        // test that PointIndex is consistent with actual point passed
        //
        for(i = 0; i <= m-1; i++)
        {
            nlserrors = nlserrors||ap::fp_neq(xy(state.pointindex,i),state.x(i));
        }
        
        //
        // calculate
        //
        if( state.needf )
        {
            v = ap::vdotproduct(&state.x(0), 1, &state.c(0), 1, ap::vlen(0,m-1));
            state.f = v;
            continue;
        }
        if( state.needfg )
        {
            v = ap::vdotproduct(&state.x(0), 1, &state.c(0), 1, ap::vlen(0,m-1));
            state.f = v;
            ap::vmove(&state.g(0), 1, &state.x(0), 1, ap::vlen(0,m-1));
            continue;
        }
        if( state.needfgh )
        {
            v = ap::vdotproduct(&state.x(0), 1, &state.c(0), 1, ap::vlen(0,m-1));
            state.f = v;
            ap::vmove(&state.g(0), 1, &state.x(0), 1, ap::vlen(0,m-1));
            for(i = 0; i <= m-1; i++)
            {
                for(j = 0; j <= m-1; j++)
                {
                    state.h(i,j) = 0;
                }
            }
            continue;
        }
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool llstestunit_test_silent()
{
    bool result;

    result = testlls(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool llstestunit_test()
{
    bool result;

    result = testlls(false);
    return result;
}





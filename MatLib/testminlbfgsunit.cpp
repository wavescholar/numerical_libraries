
 
#include <stdio.h>
#include "testminlbfgsunit.h"

static void testfunc1(minlbfgsstate& state);
static void testfunc2(minlbfgsstate& state);
static void testfunc3(minlbfgsstate& state);

bool testminlbfgs(bool silent)
{
    bool result;
    bool waserrors;
    bool referror;
    bool nonconverror;
    bool eqerror;
    bool converror;
    bool crashtest;
    bool othererrors;
    int n;
    int m;
    ap::real_1d_array x;
    ap::real_1d_array xe;
    ap::real_1d_array b;
    ap::real_1d_array xlast;
    int i;
    int j;
    double v;
    ap::real_2d_array a;
    int maxits;
    minlbfgsstate state;
    minlbfgsreport rep;
    double fprev;
    double xprev;
    double stpmax;

    waserrors = false;
    
    //
    // Reference problem
    //
    x.setbounds(0, 2);
    n = 3;
    m = 2;
    x(0) = 100*ap::randomreal()-50;
    x(1) = 100*ap::randomreal()-50;
    x(2) = 100*ap::randomreal()-50;
    minlbfgscreate(n, m, x, state);
    minlbfgssetcond(state, double(0), double(0), double(0), 0);
    while(minlbfgsiteration(state))
    {
        state.f = ap::sqr(state.x(0)-2)+ap::sqr(state.x(1))+ap::sqr(state.x(2)-state.x(0));
        state.g(0) = 2*(state.x(0)-2)+2*(state.x(0)-state.x(2));
        state.g(1) = 2*state.x(1);
        state.g(2) = 2*(state.x(2)-state.x(0));
    }
    minlbfgsresults(state, x, rep);
    referror = rep.terminationtype<=0||ap::fp_greater(fabs(x(0)-2),0.001)||ap::fp_greater(fabs(x(1)),0.001)||ap::fp_greater(fabs(x(2)-2),0.001);
    
    //
    // nonconvex problems with hard relief: we start from point with very small
    // gradient, but we need ever smaller gradient in the next step due to
    // Wolfe conditions.
    //
    nonconverror = false;
    x.setlength(1);
    n = 1;
    m = 1;
    v = -100;
    while(ap::fp_less(v,0.1))
    {
        x(0) = v;
        minlbfgscreate(n, m, x, state);
        minlbfgssetcond(state, 1.0E-9, double(0), double(0), 0);
        while(minlbfgsiteration(state))
        {
            state.f = ap::sqr(state.x(0))/(1+ap::sqr(state.x(0)));
            state.g(0) = (2*state.x(0)*(1+ap::sqr(state.x(0)))-ap::sqr(state.x(0))*2*state.x(0))/ap::sqr(1+ap::sqr(state.x(0)));
        }
        minlbfgsresults(state, x, rep);
        nonconverror = nonconverror||rep.terminationtype<=0||ap::fp_greater(fabs(x(0)),0.001);
        v = v+0.1;
    }
    
    //
    // Linear equations
    //
    eqerror = false;
    for(n = 1; n <= 10; n++)
    {
        
        //
        // Prepare task
        //
        a.setbounds(0, n-1, 0, n-1);
        x.setbounds(0, n-1);
        xe.setbounds(0, n-1);
        b.setbounds(0, n-1);
        for(i = 0; i <= n-1; i++)
        {
            xe(i) = 2*ap::randomreal()-1;
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                a(i,j) = 2*ap::randomreal()-1;
            }
            a(i,i) = a(i,i)+3*ap::sign(a(i,i));
        }
        for(i = 0; i <= n-1; i++)
        {
            v = ap::vdotproduct(&a(i, 0), 1, &xe(0), 1, ap::vlen(0,n-1));
            b(i) = v;
        }
        
        //
        // Test different M
        //
        for(m = 1; m <= n; m++)
        {
            
            //
            // Solve task
            //
            for(i = 0; i <= n-1; i++)
            {
                x(i) = 2*ap::randomreal()-1;
            }
            minlbfgscreate(n, m, x, state);
            minlbfgssetcond(state, double(0), double(0), double(0), 0);
            while(minlbfgsiteration(state))
            {
                state.f = 0;
                for(i = 0; i <= n-1; i++)
                {
                    state.g(i) = 0;
                }
                for(i = 0; i <= n-1; i++)
                {
                    v = ap::vdotproduct(&a(i, 0), 1, &state.x(0), 1, ap::vlen(0,n-1));
                    state.f = state.f+ap::sqr(v-b(i));
                    for(j = 0; j <= n-1; j++)
                    {
                        state.g(j) = state.g(j)+2*(v-b(i))*a(i,j);
                    }
                }
            }
            minlbfgsresults(state, x, rep);
            eqerror = eqerror||rep.terminationtype<=0;
            for(i = 0; i <= n-1; i++)
            {
                eqerror = eqerror||ap::fp_greater(fabs(x(i)-xe(i)),0.001);
            }
        }
    }
    
    //
    // Testing convergence properties
    //
    converror = false;
    x.setbounds(0, 2);
    n = 3;
    m = 2;
    for(i = 0; i <= 2; i++)
    {
        x(i) = 6*ap::randomreal()-3;
    }
    minlbfgscreate(n, m, x, state);
    minlbfgssetcond(state, 0.001, double(0), double(0), 0);
    while(minlbfgsiteration(state))
    {
        testfunc3(state);
    }
    minlbfgsresults(state, x, rep);
    converror = converror||rep.terminationtype!=4;
    for(i = 0; i <= 2; i++)
    {
        x(i) = 6*ap::randomreal()-3;
    }
    minlbfgscreate(n, m, x, state);
    minlbfgssetcond(state, double(0), 0.001, double(0), 0);
    while(minlbfgsiteration(state))
    {
        testfunc3(state);
    }
    minlbfgsresults(state, x, rep);
    converror = converror||rep.terminationtype!=1;
    for(i = 0; i <= 2; i++)
    {
        x(i) = 6*ap::randomreal()-3;
    }
    minlbfgscreate(n, m, x, state);
    minlbfgssetcond(state, double(0), double(0), 0.001, 0);
    while(minlbfgsiteration(state))
    {
        testfunc3(state);
    }
    minlbfgsresults(state, x, rep);
    converror = converror||rep.terminationtype!=2;
    for(i = 0; i <= 2; i++)
    {
        x(i) = 2*ap::randomreal()-1;
    }
    minlbfgscreate(n, m, x, state);
    minlbfgssetcond(state, double(0), double(0), double(0), 10);
    while(minlbfgsiteration(state))
    {
        testfunc3(state);
    }
    minlbfgsresults(state, x, rep);
    converror = converror||rep.terminationtype!=5||rep.iterationscount!=10;
    
    //
    // Crash test: too many iterations on a simple tasks
    // May fail when encounter zero step, underflow or something like that
    //
    crashtest = false;
    x.setbounds(0, 2);
    n = 3;
    m = 2;
    maxits = 10000;
    for(i = 0; i <= 2; i++)
    {
        x(i) = 6*ap::randomreal()-3;
    }
    minlbfgscreate(n, m, x, state);
    minlbfgssetcond(state, double(0), double(0), double(0), maxits);
    while(minlbfgsiteration(state))
    {
        state.f = ap::sqr(exp(state.x(0))-2)+ap::sqr(state.x(1))+ap::sqr(state.x(2)-state.x(0));
        state.g(0) = 2*(exp(state.x(0))-2)*exp(state.x(0))+2*(state.x(0)-state.x(2));
        state.g(1) = 2*state.x(1);
        state.g(2) = 2*(state.x(2)-state.x(0));
    }
    minlbfgsresults(state, x, rep);
    crashtest = crashtest||rep.terminationtype<=0;
    
    //
    // Other properties:
    // 1. test reports (F should form monotone sequence)
    // 2. test maximum step
    //
    othererrors = false;
    n = 50;
    m = 2;
    x.setlength(n);
    xlast.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        x(i) = 1;
    }
    minlbfgscreate(n, m, x, state);
    minlbfgssetcond(state, double(0), double(0), double(0), 100);
    minlbfgssetxrep(state, true);
    fprev = ap::maxrealnumber;
    while(minlbfgsiteration(state))
    {
        if( state.needfg )
        {
            state.f = 0;
            for(i = 0; i <= n-1; i++)
            {
                state.f = state.f+ap::sqr((1+i)*state.x(i));
                state.g(i) = 2*(1+i)*state.x(i);
            }
        }
        if( state.xupdated )
        {
            othererrors = othererrors||ap::fp_greater(state.f,fprev);
            if( ap::fp_eq(fprev,ap::maxrealnumber) )
            {
                for(i = 0; i <= n-1; i++)
                {
                    othererrors = othererrors||ap::fp_neq(state.x(i),x(i));
                }
            }
            fprev = state.f;
            ap::vmove(&xlast(0), 1, &state.x(0), 1, ap::vlen(0,n-1));
        }
    }
    minlbfgsresults(state, x, rep);
    for(i = 0; i <= n-1; i++)
    {
        othererrors = othererrors||ap::fp_neq(x(i),xlast(i));
    }
    n = 1;
    m = 1;
    x.setlength(n);
    x(0) = 100;
    stpmax = 0.05+0.05*ap::randomreal();
    minlbfgscreate(n, m, x, state);
    minlbfgssetcond(state, 1.0E-9, double(0), double(0), 0);
    minlbfgssetstpmax(state, stpmax);
    minlbfgssetxrep(state, true);
    xprev = x(0);
    while(minlbfgsiteration(state))
    {
        if( state.needfg )
        {
            state.f = exp(state.x(0))+exp(-state.x(0));
            state.g(0) = exp(state.x(0))-exp(-state.x(0));
            othererrors = othererrors||ap::fp_greater(fabs(state.x(0)-xprev),(1+sqrt(ap::machineepsilon))*stpmax);
        }
        if( state.xupdated )
        {
            othererrors = othererrors||ap::fp_greater(fabs(state.x(0)-xprev),(1+sqrt(ap::machineepsilon))*stpmax);
            xprev = state.x(0);
        }
    }
    
    //
    // end
    //
    waserrors = referror||nonconverror||eqerror||converror||crashtest||othererrors;
    if( !silent )
    {
        printf("TESTING L-BFGS OPTIMIZATION\n");
        printf("REFERENCE PROBLEM:                        ");
        if( referror )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("NON-CONVEX PROBLEM:                       ");
        if( nonconverror )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("LINEAR EQUATIONS:                         ");
        if( eqerror )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("CONVERGENCE PROPERTIES:                   ");
        if( converror )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("CRASH TEST:                               ");
        if( crashtest )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("OTHER PROPERTIES:                         ");
        if( othererrors )
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
    result = !waserrors;
    return result;
}


/*************************************************************************
Calculate test function #1

It may show very interesting behavior when optimized with 'x[0]>=ln(2)'
constraint.
*************************************************************************/
static void testfunc1(minlbfgsstate& state)
{

    if( ap::fp_less(state.x(0),100) )
    {
        state.f = ap::sqr(exp(state.x(0))-2)+ap::sqr(state.x(1))+ap::sqr(state.x(2)-state.x(0));
        state.g(0) = 2*(exp(state.x(0))-2)*exp(state.x(0))+2*(state.x(0)-state.x(2));
        state.g(1) = 2*state.x(1);
        state.g(2) = 2*(state.x(2)-state.x(0));
    }
    else
    {
        state.f = sqrt(ap::maxrealnumber);
        state.g(0) = sqrt(ap::maxrealnumber);
        state.g(1) = 0;
        state.g(2) = 0;
    }
}


/*************************************************************************
Calculate test function #2

Simple variation of #1, much more nonlinear, which makes unlikely premature
convergence of algorithm .
*************************************************************************/
static void testfunc2(minlbfgsstate& state)
{

    if( ap::fp_less(state.x(0),100) )
    {
        state.f = ap::sqr(exp(state.x(0))-2)+ap::sqr(ap::sqr(state.x(1)))+ap::sqr(state.x(2)-state.x(0));
        state.g(0) = 2*(exp(state.x(0))-2)*exp(state.x(0))+2*(state.x(0)-state.x(2));
        state.g(1) = 4*state.x(1)*ap::sqr(state.x(1));
        state.g(2) = 2*(state.x(2)-state.x(0));
    }
    else
    {
        state.f = sqrt(ap::maxrealnumber);
        state.g(0) = sqrt(ap::maxrealnumber);
        state.g(1) = 0;
        state.g(2) = 0;
    }
}


/*************************************************************************
Calculate test function #3

Simple variation of #1, much more nonlinear, with non-zero value at minimum.
It achieve two goals:
* makes unlikely premature convergence of algorithm .
* solves some issues with EpsF stopping condition which arise when
  F(minimum) is zero

*************************************************************************/
static void testfunc3(minlbfgsstate& state)
{
    double s;

    s = 0.001;
    if( ap::fp_less(state.x(0),100) )
    {
        state.f = ap::sqr(exp(state.x(0))-2)+ap::sqr(ap::sqr(state.x(1))+s)+ap::sqr(state.x(2)-state.x(0));
        state.g(0) = 2*(exp(state.x(0))-2)*exp(state.x(0))+2*(state.x(0)-state.x(2));
        state.g(1) = 2*(ap::sqr(state.x(1))+s)*2*state.x(1);
        state.g(2) = 2*(state.x(2)-state.x(0));
    }
    else
    {
        state.f = sqrt(ap::maxrealnumber);
        state.g(0) = sqrt(ap::maxrealnumber);
        state.g(1) = 0;
        state.g(2) = 0;
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testminlbfgsunit_test_silent()
{
    bool result;

    result = testminlbfgs(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testminlbfgsunit_test()
{
    bool result;

    result = testminlbfgs(false);
    return result;
}





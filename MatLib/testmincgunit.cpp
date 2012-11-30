
 
#include <stdio.h>
#include "testmincgunit.h"

static void testfunc1(mincgstate& state);
static void testfunc2(mincgstate& state);
static void testfunc3(mincgstate& state);

bool testmincg(bool silent)
{
    bool result;
    bool waserrors;
    bool referror;
    bool eqerror;
    bool linerror1;
    bool linerror2;
    bool converror;
    bool othererrors;
    int n;
    ap::real_1d_array x;
    ap::real_1d_array xe;
    ap::real_1d_array b;
    ap::real_1d_array xlast;
    double fprev;
    double xprev;
    double stpmax;
    int i;
    int j;
    double v;
    ap::real_2d_array a;
    mincgstate state;
    mincgreport rep;
    int cgtype;

    waserrors = false;
    referror = false;
    linerror1 = false;
    linerror2 = false;
    eqerror = false;
    converror = false;
    othererrors = false;
    for(cgtype = 0; cgtype <= 1; cgtype++)
    {
        
        //
        // Reference problem
        //
        x.setbounds(0, 2);
        n = 3;
        x(0) = 100*ap::randomreal()-50;
        x(1) = 100*ap::randomreal()-50;
        x(2) = 100*ap::randomreal()-50;
        mincgcreate(n, x, state);
        mincgsetcgtype(state, cgtype);
        while(mincgiteration(state))
        {
            state.f = ap::sqr(state.x(0)-2)+ap::sqr(state.x(1))+ap::sqr(state.x(2)-state.x(0));
            state.g(0) = 2*(state.x(0)-2)+2*(state.x(0)-state.x(2));
            state.g(1) = 2*state.x(1);
            state.g(2) = 2*(state.x(2)-state.x(0));
        }
        mincgresults(state, x, rep);
        referror = referror||rep.terminationtype<=0||ap::fp_greater(fabs(x(0)-2),0.001)||ap::fp_greater(fabs(x(1)),0.001)||ap::fp_greater(fabs(x(2)-2),0.001);
        
        //
        // 1D problem #1
        //
        x.setbounds(0, 0);
        n = 1;
        x(0) = 100*ap::randomreal()-50;
        mincgcreate(n, x, state);
        mincgsetcgtype(state, cgtype);
        while(mincgiteration(state))
        {
            state.f = -cos(state.x(0));
            state.g(0) = sin(state.x(0));
        }
        mincgresults(state, x, rep);
        linerror1 = linerror1||rep.terminationtype<=0||ap::fp_greater(fabs(x(0)/ap::pi()-ap::round(x(0)/ap::pi())),0.001);
        
        //
        // 1D problem #2
        //
        x.setbounds(0, 0);
        n = 1;
        x(0) = 100*ap::randomreal()-50;
        mincgcreate(n, x, state);
        mincgsetcgtype(state, cgtype);
        while(mincgiteration(state))
        {
            state.f = ap::sqr(state.x(0))/(1+ap::sqr(state.x(0)));
            state.g(0) = (2*state.x(0)*(1+ap::sqr(state.x(0)))-ap::sqr(state.x(0))*2*state.x(0))/ap::sqr(1+ap::sqr(state.x(0)));
        }
        mincgresults(state, x, rep);
        linerror2 = linerror2||rep.terminationtype<=0||ap::fp_greater(fabs(x(0)),0.001);
        
        //
        // Linear equations
        //
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
            // Solve task
            //
            for(i = 0; i <= n-1; i++)
            {
                x(i) = 2*ap::randomreal()-1;
            }
            mincgcreate(n, x, state);
            mincgsetcgtype(state, cgtype);
            while(mincgiteration(state))
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
            mincgresults(state, x, rep);
            eqerror = eqerror||rep.terminationtype<=0;
            for(i = 0; i <= n-1; i++)
            {
                eqerror = eqerror||ap::fp_greater(fabs(x(i)-xe(i)),0.001);
            }
        }
        
        //
        // Testing convergence properties
        //
        x.setbounds(0, 2);
        n = 3;
        for(i = 0; i <= 2; i++)
        {
            x(i) = 6*ap::randomreal()-3;
        }
        mincgcreate(n, x, state);
        mincgsetcond(state, 0.001, 0.0, 0.0, 0);
        mincgsetcgtype(state, cgtype);
        while(mincgiteration(state))
        {
            testfunc3(state);
        }
        mincgresults(state, x, rep);
        converror = converror||rep.terminationtype!=4;
        for(i = 0; i <= 2; i++)
        {
            x(i) = 6*ap::randomreal()-3;
        }
        mincgcreate(n, x, state);
        mincgsetcond(state, 0.0, 0.001, 0.0, 0);
        mincgsetcgtype(state, cgtype);
        while(mincgiteration(state))
        {
            testfunc3(state);
        }
        mincgresults(state, x, rep);
        converror = converror||rep.terminationtype!=1;
        for(i = 0; i <= 2; i++)
        {
            x(i) = 6*ap::randomreal()-3;
        }
        mincgcreate(n, x, state);
        mincgsetcond(state, 0.0, 0.0, 0.001, 0);
        mincgsetcgtype(state, cgtype);
        while(mincgiteration(state))
        {
            testfunc3(state);
        }
        mincgresults(state, x, rep);
        converror = converror||rep.terminationtype!=2;
        for(i = 0; i <= 2; i++)
        {
            x(i) = 2*ap::randomreal()-1;
        }
        mincgcreate(n, x, state);
        mincgsetcond(state, 0.0, 0.0, 0.0, 10);
        mincgsetcgtype(state, cgtype);
        while(mincgiteration(state))
        {
            testfunc3(state);
        }
        mincgresults(state, x, rep);
        converror = converror||!(rep.terminationtype==5&&rep.iterationscount==10||rep.terminationtype==7);
        
        //
        // Other properties:
        // 1. test reports (F should form monotone sequence)
        // 2. test maximum step
        //
        n = 50;
        x.setlength(n);
        xlast.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            x(i) = 1;
        }
        mincgcreate(n, x, state);
        mincgsetcond(state, double(0), double(0), double(0), 100);
        mincgsetxrep(state, true);
        fprev = ap::maxrealnumber;
        while(mincgiteration(state))
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
        mincgresults(state, x, rep);
        for(i = 0; i <= n-1; i++)
        {
            othererrors = othererrors||ap::fp_neq(x(i),xlast(i));
        }
        n = 1;
        x.setlength(n);
        x(0) = 100;
        stpmax = 0.05+0.05*ap::randomreal();
        mincgcreate(n, x, state);
        mincgsetcond(state, 1.0E-9, double(0), double(0), 0);
        mincgsetstpmax(state, stpmax);
        mincgsetxrep(state, true);
        xprev = x(0);
        while(mincgiteration(state))
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
    }
    
    //
    // end
    //
    waserrors = referror||eqerror||linerror1||linerror2||converror||othererrors;
    if( !silent )
    {
        printf("TESTING CG OPTIMIZATION\n");
        printf("REFERENCE PROBLEM:                        ");
        if( referror )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("LIN-1 PROBLEM:                            ");
        if( linerror1 )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("LIN-2 PROBLEM:                            ");
        if( linerror2 )
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
*************************************************************************/
static void testfunc1(mincgstate& state)
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
static void testfunc2(mincgstate& state)
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
static void testfunc3(mincgstate& state)
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
bool testmincgunit_test_silent()
{
    bool result;

    result = testmincg(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testmincgunit_test()
{
    bool result;

    result = testmincg(false);
    return result;
}





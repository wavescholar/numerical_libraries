
 
#include <stdio.h>
#include "testodesolverunit.h"

static void unset2d(ap::real_2d_array& x);
static void unset1d(ap::real_1d_array& x);
static void unsetrep(odesolverreport& rep);

/*************************************************************************
Test
*************************************************************************/
bool testodesolver(bool silent)
{
    bool result;
    int passcount;
    bool curerrors;
    bool rkckerrors;
    bool waserrors;
    ap::real_1d_array xtbl;
    ap::real_2d_array ytbl;
    odesolverreport rep;
    ap::real_1d_array xg;
    ap::real_1d_array y;
    double h;
    double eps;
    int solver;
    int pass;
    int mynfev;
    double v;
    int n;
    int m;
    int m2;
    int i;
    double err;
    odesolverstate state;

    rkckerrors = false;
    waserrors = false;
    passcount = 10;
    
    //
    // simple test: just A*sin(x)+B*cos(x)
    //
    ap::ap_error::make_assertion(passcount>=2, "");
    for(pass = 0; pass <= passcount-1; pass++)
    {
        for(solver = 0; solver <= 0; solver++)
        {
            
            //
            // prepare
            //
            h = 1.0E-2;
            eps = 1.0E-5;
            if( pass%2==0 )
            {
                eps = -eps;
            }
            y.setlength(2);
            for(i = 0; i <= 1; i++)
            {
                y(i) = 2*ap::randomreal()-1;
            }
            m = 2+ap::randominteger(10);
            xg.setlength(m);
            xg(0) = (m-1)*ap::randomreal();
            for(i = 1; i <= m-1; i++)
            {
                xg(i) = xg(i-1)+ap::randomreal();
            }
            v = 2*ap::pi()/(xg(m-1)-xg(0));
            ap::vmul(&xg(0), 1, ap::vlen(0,m-1), v);
            if( ap::fp_greater(ap::randomreal(),0.5) )
            {
                ap::vmul(&xg(0), 1, ap::vlen(0,m-1), -1);
            }
            mynfev = 0;
            
            //
            // choose solver
            //
            if( solver==0 )
            {
                odesolverrkck(y, 2, xg, m, eps, h, state);
            }
            
            //
            // solve
            //
            while(odesolveriteration(state))
            {
                state.dy(0) = state.y(1);
                state.dy(1) = -state.y(0);
                mynfev = mynfev+1;
            }
            odesolverresults(state, m2, xtbl, ytbl, rep);
            
            //
            // check results
            //
            curerrors = false;
            if( rep.terminationtype<=0 )
            {
                curerrors = true;
            }
            else
            {
                curerrors = curerrors||m2!=m;
                err = 0;
                for(i = 0; i <= m-1; i++)
                {
                    err = ap::maxreal(err, fabs(ytbl(i,0)-(+y(0)*cos(xtbl(i)-xtbl(0))+y(1)*sin(xtbl(i)-xtbl(0)))));
                    err = ap::maxreal(err, fabs(ytbl(i,1)-(-y(0)*sin(xtbl(i)-xtbl(0))+y(1)*cos(xtbl(i)-xtbl(0)))));
                }
                curerrors = curerrors||ap::fp_greater(err,10*fabs(eps));
                curerrors = curerrors||mynfev!=rep.nfev;
            }
            if( solver==0 )
            {
                rkckerrors = rkckerrors||curerrors;
            }
        }
    }
    
    //
    // another test:
    //
    //     y(0)   = 0
    //     dy/dx  = f(x,y)
    //     f(x,y) = 0,   x<1
    //              x-1, x>=1
    //
    // with BOTH absolute and fractional tolerances.
    // Starting from zero will be real challenge for
    // fractional tolerance.
    //
    ap::ap_error::make_assertion(passcount>=2, "");
    for(pass = 0; pass <= passcount-1; pass++)
    {
        h = 1.0E-4;
        eps = 1.0E-4;
        if( pass%2==0 )
        {
            eps = -eps;
        }
        y.setlength(1);
        y(0) = 0;
        m = 21;
        xg.setlength(m);
        for(i = 0; i <= m-1; i++)
        {
            xg(i) = double(2*i)/double(m-1);
        }
        mynfev = 0;
        odesolverrkck(y, 1, xg, m, eps, h, state);
        while(odesolveriteration(state))
        {
            state.dy(0) = ap::maxreal(state.x-1, double(0));
            mynfev = mynfev+1;
        }
        odesolverresults(state, m2, xtbl, ytbl, rep);
        if( rep.terminationtype<=0 )
        {
            rkckerrors = true;
        }
        else
        {
            rkckerrors = rkckerrors||m2!=m;
            err = 0;
            for(i = 0; i <= m-1; i++)
            {
                err = ap::maxreal(err, fabs(ytbl(i,0)-ap::sqr(ap::maxreal(xg(i)-1, double(0)))/2));
            }
            rkckerrors = rkckerrors||ap::fp_greater(err,fabs(eps));
            rkckerrors = rkckerrors||mynfev!=rep.nfev;
        }
    }
    
    //
    // end
    //
    waserrors = rkckerrors;
    if( !silent )
    {
        printf("TESTING ODE SOLVER\n");
        printf("* RK CASH-KARP:                           ");
        if( rkckerrors )
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
Unsets real matrix
*************************************************************************/
static void unset2d(ap::real_2d_array& x)
{

    x.setlength(1, 1);
    x(0,0) = 2*ap::randomreal()-1;
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
Unsets report
*************************************************************************/
static void unsetrep(odesolverreport& rep)
{

    rep.nfev = 0;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testodesolverunit_test_silent()
{
    bool result;

    result = testodesolver(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testodesolverunit_test()
{
    bool result;

    result = testodesolver(false);
    return result;
}





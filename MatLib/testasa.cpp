
 
#include <stdio.h>
#include "testasa.h"

static void testfunc1(minasastate& state);
static void testfunc2(minasastate& state);
static void testfunc3(minasastate& state);
static void checkbounds(const ap::real_1d_array& x,
     const ap::real_1d_array& bndl,
     const ap::real_1d_array& bndu,
     int n,
     bool& err);
static double asaboundval(double x, double b1, double b2);

bool testminasa(bool silent)
{
    bool result;
    bool waserrors;
    bool referror;
    bool converror;
    bool othererrors;
    int n;
    ap::real_1d_array x;
    ap::real_1d_array xe;
    ap::real_1d_array c;
    ap::real_1d_array bndl;
    ap::real_1d_array bndu;
    ap::real_1d_array xlast;
    double fprev;
    double xprev;
    double stpmax;
    int i;
    int j;
    double v;
    double s;
    double tol;
    int algotype;
    ap::real_2d_array a;
    minasastate state;
    minasareport rep;

    waserrors = false;
    referror = false;
    converror = false;
    othererrors = false;
    
    //
    // Different algorithms
    //
    for(algotype = -1; algotype <= 1; algotype++)
    {
        
        //
        // reference problem, simple convex optimization
        //
        for(n = 1; n <= 5; n++)
        {
            
            //
            // min(x'*diag(c)*x) on a random box
            //
            x.setlength(n);
            xe.setlength(n);
            c.setlength(n);
            bndl.setlength(n);
            bndu.setlength(n);
            for(i = 0; i <= n-1; i++)
            {
                c(i) = 1+ap::randomreal();
                xe(i) = 4*ap::randomreal()-2;
                bndl(i) = -ap::maxreal(ap::randomreal(), 0.2);
                bndu(i) = +ap::maxreal(ap::randomreal(), 0.2);
                x(i) = 0.5*(bndl(i)+bndu(i));
            }
            tol = 0.001;
            minasacreate(n, x, bndl, bndu, state);
            minasasetcond(state, tol, 0.0, 0.0, 0);
            minasasetalgorithm(state, algotype);
            while(minasaiteration(state))
            {
                checkbounds(state.x, bndl, bndu, n, othererrors);
                state.f = 0;
                for(i = 0; i <= n-1; i++)
                {
                    state.f = state.f+c(i)*ap::sqr(state.x(i)-xe(i));
                    state.g(i) = 2*c(i)*(state.x(i)-xe(i));
                }
            }
            minasaresults(state, x, rep);
            referror = referror||rep.terminationtype<=0;
            for(i = 0; i <= n-1; i++)
            {
                referror = referror||ap::fp_greater(fabs(asaboundval(xe(i), bndl(i), bndu(i))-x(i)),0.01);
            }
        }
        
        //
        // reference problem 2: non-convex optimization on [-2,2] x [1,2]
        //
        // A saddle function is minimized:
        // * stationary point [0,0] (non-feasible)
        // * constrained minimum [-2,2].
        // * starting point [+2,2]
        //
        // Path from start to end may be very complex, with multiple changes
        // in active constraints, so it is interesting task for our method.
        //
        // Scale parameter is used to make optimization more interesting
        // during GPA runs.
        //
        x.setlength(2);
        bndl.setlength(2);
        bndu.setlength(2);
        bndl(0) = -2;
        bndu(0) = 2;
        x(0) = 2;
        bndl(1) = 1;
        bndu(1) = 2;
        x(1) = 2;
        tol = 0.001;
        s = 0.01;
        minasacreate(2, x, bndl, bndu, state);
        minasasetcond(state, tol, 0.0, 0.0, 0);
        minasasetalgorithm(state, algotype);
        while(minasaiteration(state))
        {
            checkbounds(state.x, bndl, bndu, 2, othererrors);
            state.f = s*(ap::sqr(state.x(0)+state.x(1))-ap::sqr(state.x(0)-state.x(1)));
            state.g(0) = s*(2*(state.x(0)+state.x(1))-2*(state.x(0)-state.x(1)));
            state.g(1) = s*(2*(state.x(0)+state.x(1))+2*(state.x(0)-state.x(1)));
        }
        minasaresults(state, x, rep);
        referror = referror||rep.terminationtype<=0||ap::fp_greater(fabs(state.x(0)+2),0.01)||ap::fp_greater(fabs(state.x(1)-2),0.01);
        
        //
        // function #1 with 'x[0]>=ln(2)' constraint.
        // may show very interesting behavior.
        //
        x.setlength(3);
        bndl.setlength(3);
        bndu.setlength(3);
        n = 3;
        for(i = 0; i <= 2; i++)
        {
            bndl(i) = -10000;
            bndu(i) = +10000;
        }
        bndl(0) = log(double(2));
        for(i = 0; i <= 2; i++)
        {
            x(i) = 3*ap::randomreal()+3;
        }
        minasacreate(n, x, bndl, bndu, state);
        minasasetcond(state, 0.0000001, 0.0, 0.0, 0);
        minasasetalgorithm(state, algotype);
        while(minasaiteration(state))
        {
            checkbounds(state.x, bndl, bndu, n, othererrors);
            testfunc1(state);
        }
        minasaresults(state, x, rep);
        referror = referror||rep.terminationtype<=0;
        referror = referror||ap::fp_greater(fabs(x(0)-log(double(2))),0.05);
        referror = referror||ap::fp_greater(fabs(x(1)),0.05);
        referror = referror||ap::fp_greater(fabs(x(2)-log(double(2))),0.05);
        
        //
        // Testing convergence properties
        //
        x.setlength(3);
        bndl.setlength(3);
        bndu.setlength(3);
        n = 3;
        for(i = 0; i <= 2; i++)
        {
            bndl(i) = -10000;
            bndu(i) = +10000;
        }
        bndl(0) = log(double(2));
        for(i = 0; i <= 2; i++)
        {
            x(i) = 3*ap::randomreal()+3;
        }
        minasacreate(n, x, bndl, bndu, state);
        minasasetcond(state, 0.001, 0.0, 0.0, 0);
        minasasetalgorithm(state, algotype);
        while(minasaiteration(state))
        {
            checkbounds(state.x, bndl, bndu, n, othererrors);
            testfunc3(state);
        }
        minasaresults(state, x, rep);
        converror = converror||rep.terminationtype!=4;
        for(i = 0; i <= 2; i++)
        {
            x(i) = 3*ap::randomreal()+3;
        }
        minasacreate(n, x, bndl, bndu, state);
        minasasetcond(state, 0.0, 0.001, 0.0, 0);
        minasasetalgorithm(state, algotype);
        while(minasaiteration(state))
        {
            checkbounds(state.x, bndl, bndu, n, othererrors);
            testfunc3(state);
        }
        minasaresults(state, x, rep);
        converror = converror||rep.terminationtype!=1;
        for(i = 0; i <= 2; i++)
        {
            x(i) = 3*ap::randomreal()+3;
        }
        minasacreate(n, x, bndl, bndu, state);
        minasasetcond(state, 0.0, 0.0, 0.001, 0);
        minasasetalgorithm(state, algotype);
        while(minasaiteration(state))
        {
            checkbounds(state.x, bndl, bndu, n, othererrors);
            testfunc3(state);
        }
        minasaresults(state, x, rep);
        converror = converror||rep.terminationtype!=2;
        for(i = 0; i <= 2; i++)
        {
            x(i) = 3*ap::randomreal()+3;
        }
        minasacreate(n, x, bndl, bndu, state);
        minasasetcond(state, 0.0, 0.0, 0.0, 3);
        minasasetalgorithm(state, algotype);
        while(minasaiteration(state))
        {
            checkbounds(state.x, bndl, bndu, n, othererrors);
            testfunc3(state);
        }
        minasaresults(state, x, rep);
        converror = converror||!(rep.terminationtype==5&&rep.iterationscount==3||rep.terminationtype==7);
        
        //
        // Other properties
        //
        //
        // Other properties:
        // 1. test reports (F should form monotone sequence)
        // 2. test maximum step
        //
        n = 50;
        x.setlength(n);
        xlast.setlength(n);
        bndl.setlength(n);
        bndu.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            x(i) = 1;
            xlast(i) = ap::randomreal();
            bndl(i) = -100000;
            bndu(i) = +100000;
        }
        minasacreate(n, x, bndl, bndu, state);
        minasasetcond(state, double(0), double(0), double(0), 100);
        minasasetxrep(state, true);
        fprev = ap::maxrealnumber;
        while(minasaiteration(state))
        {
            checkbounds(state.x, bndl, bndu, n, othererrors);
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
        minasaresults(state, x, rep);
        for(i = 0; i <= n-1; i++)
        {
            othererrors = othererrors||ap::fp_neq(x(i),xlast(i));
        }
        n = 1;
        x.setlength(n);
        bndl.setlength(n);
        bndu.setlength(n);
        x(0) = 100;
        bndl(0) = -1000000;
        bndu(0) = +1000000;
        stpmax = 0.05+0.05*ap::randomreal();
        minasacreate(n, x, bndl, bndu, state);
        minasasetcond(state, 1.0E-9, double(0), double(0), 0);
        minasasetstpmax(state, stpmax);
        minasasetxrep(state, true);
        xprev = x(0);
        while(minasaiteration(state))
        {
            checkbounds(state.x, bndl, bndu, n, othererrors);
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
    waserrors = referror||converror||othererrors;
    if( !silent )
    {
        printf("TESTING ASA OPTIMIZATION\n");
        printf("REFERENCE PROBLEMS:                       ");
        if( referror )
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

It may show very interesting behavior when optimized with 'x[0]>=ln(2)'
constraint.
*************************************************************************/
static void testfunc1(minasastate& state)
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
static void testfunc2(minasastate& state)
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
static void testfunc3(minasastate& state)
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
Checks that X is bounded with respect to BndL/BndU.

If it is not, True is assigned to the Err variable (which is not changed
otherwise).
*************************************************************************/
static void checkbounds(const ap::real_1d_array& x,
     const ap::real_1d_array& bndl,
     const ap::real_1d_array& bndu,
     int n,
     bool& err)
{
    int i;

    for(i = 0; i <= n-1; i++)
    {
        if( ap::fp_less(x(i),bndl(i))||ap::fp_greater(x(i),bndu(i)) )
        {
            err = true;
        }
    }
}


/*************************************************************************
'bound' value: map X to [B1,B2]
*************************************************************************/
static double asaboundval(double x, double b1, double b2)
{
    double result;

    if( ap::fp_less_eq(x,b1) )
    {
        result = b1;
        return result;
    }
    if( ap::fp_greater_eq(x,b2) )
    {
        result = b2;
        return result;
    }
    result = x;
    return result;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testasa_test_silent()
{
    bool result;

    result = testminasa(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testasa_test()
{
    bool result;

    result = testminasa(false);
    return result;
}





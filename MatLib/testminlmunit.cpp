
 
#include <stdio.h>
#include "testminlmunit.h"

static bool rkindvsstatecheck(int rkind, const minlmstate& state);

bool testminlm(bool silent)
{
    bool result;
    bool waserrors;
    bool referror;
    bool lin1error;
    bool lin2error;
    bool eqerror;
    bool converror;
    bool scerror;
    bool othererrors;
    int rkind;
    int ckind;
    double epsf;
    double epsx;
    double epsg;
    int maxits;
    int n;
    int m;
    ap::real_1d_array x;
    ap::real_1d_array xe;
    ap::real_1d_array b;
    ap::real_1d_array xlast;
    int i;
    int j;
    int k;
    double v;
    double s;
    double stpmax;
    ap::real_2d_array a;
    double fprev;
    double xprev;
    minlmstate state;
    minlmreport rep;

    waserrors = false;
    referror = false;
    lin1error = false;
    lin2error = false;
    eqerror = false;
    converror = false;
    scerror = false;
    othererrors = false;
    
    //
    // Reference problem.
    // RKind is a algorithm selector:
    // * 0 = FJ
    // * 1 = FGJ
    // * 2 = FGH
    //
    x.setbounds(0, 2);
    n = 3;
    m = 3;
    for(rkind = 0; rkind <= 2; rkind++)
    {
        x(0) = 100*ap::randomreal()-50;
        x(1) = 100*ap::randomreal()-50;
        x(2) = 100*ap::randomreal()-50;
        if( rkind==0 )
        {
            minlmcreatefj(n, m, x, state);
        }
        if( rkind==1 )
        {
            minlmcreatefgj(n, m, x, state);
        }
        if( rkind==2 )
        {
            minlmcreatefgh(n, x, state);
        }
        while(minlmiteration(state))
        {
            
            //
            // (x-2)^2 + y^2 + (z-x)^2
            //
            state.f = ap::sqr(state.x(0)-2)+ap::sqr(state.x(1))+ap::sqr(state.x(2)-state.x(0));
            if( state.needfg||state.needfgh )
            {
                state.g(0) = 2*(state.x(0)-2)+2*(state.x(0)-state.x(2));
                state.g(1) = 2*state.x(1);
                state.g(2) = 2*(state.x(2)-state.x(0));
            }
            if( state.needfij )
            {
                state.fi(0) = state.x(0)-2;
                state.fi(1) = state.x(1);
                state.fi(2) = state.x(2)-state.x(0);
                state.j(0,0) = 1;
                state.j(0,1) = 0;
                state.j(0,2) = 0;
                state.j(1,0) = 0;
                state.j(1,1) = 1;
                state.j(1,2) = 0;
                state.j(2,0) = -1;
                state.j(2,1) = 0;
                state.j(2,2) = 1;
            }
            if( state.needfgh )
            {
                state.h(0,0) = 4;
                state.h(0,1) = 0;
                state.h(0,2) = -2;
                state.h(1,0) = 0;
                state.h(1,1) = 2;
                state.h(1,2) = 0;
                state.h(2,0) = -2;
                state.h(2,1) = 0;
                state.h(2,2) = 2;
            }
            scerror = scerror||!rkindvsstatecheck(rkind, state);
        }
        minlmresults(state, x, rep);
        referror = referror||rep.terminationtype<=0||ap::fp_greater(fabs(x(0)-2),0.001)||ap::fp_greater(fabs(x(1)),0.001)||ap::fp_greater(fabs(x(2)-2),0.001);
    }
    
    //
    // 1D problem #1
    //
    for(rkind = 0; rkind <= 2; rkind++)
    {
        x.setlength(1);
        n = 1;
        m = 1;
        x(0) = 100*ap::randomreal()-50;
        if( rkind==0 )
        {
            minlmcreatefj(n, m, x, state);
        }
        if( rkind==1 )
        {
            minlmcreatefgj(n, m, x, state);
        }
        if( rkind==2 )
        {
            minlmcreatefgh(n, x, state);
        }
        while(minlmiteration(state))
        {
            state.f = ap::sqr(sin(state.x(0)));
            if( state.needfg||state.needfgh )
            {
                state.g(0) = 2*sin(state.x(0))*cos(state.x(0));
            }
            if( state.needfij )
            {
                state.fi(0) = sin(state.x(0));
                state.j(0,0) = cos(state.x(0));
            }
            if( state.needfgh )
            {
                state.h(0,0) = 2*(cos(state.x(0))*cos(state.x(0))-sin(state.x(0))*sin(state.x(0)));
            }
            scerror = scerror||!rkindvsstatecheck(rkind, state);
        }
        minlmresults(state, x, rep);
        lin1error = rep.terminationtype<=0||ap::fp_greater(fabs(x(0)/ap::pi()-ap::round(x(0)/ap::pi())),0.001);
    }
    
    //
    // Linear equations
    //
    for(n = 1; n <= 10; n++)
    {
        
        //
        // Prepare task
        //
        rmatrixrndcond(n, double(100), a);
        x.setlength(n);
        xe.setlength(n);
        b.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            xe(i) = 2*ap::randomreal()-1;
        }
        for(i = 0; i <= n-1; i++)
        {
            v = ap::vdotproduct(&a(i, 0), 1, &xe(0), 1, ap::vlen(0,n-1));
            b(i) = v;
        }
        
        //
        // Test different RKind
        //
        for(rkind = 0; rkind <= 2; rkind++)
        {
            
            //
            // Solve task
            //
            for(i = 0; i <= n-1; i++)
            {
                x(i) = 2*ap::randomreal()-1;
            }
            if( rkind==0 )
            {
                minlmcreatefj(n, n, x, state);
            }
            if( rkind==1 )
            {
                minlmcreatefgj(n, n, x, state);
            }
            if( rkind==2 )
            {
                minlmcreatefgh(n, x, state);
            }
            while(minlmiteration(state))
            {
                if( state.needf||state.needfg||state.needfgh )
                {
                    state.f = 0;
                }
                if( state.needfg||state.needfgh )
                {
                    for(i = 0; i <= n-1; i++)
                    {
                        state.g(i) = 0;
                    }
                }
                if( state.needfgh )
                {
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= n-1; j++)
                        {
                            state.h(i,j) = 0;
                        }
                    }
                }
                for(i = 0; i <= n-1; i++)
                {
                    v = ap::vdotproduct(&a(i, 0), 1, &state.x(0), 1, ap::vlen(0,n-1));
                    if( state.needf||state.needfg||state.needfgh )
                    {
                        state.f = state.f+ap::sqr(v-b(i));
                    }
                    if( state.needfg||state.needfgh )
                    {
                        for(j = 0; j <= n-1; j++)
                        {
                            state.g(j) = state.g(j)+2*(v-b(i))*a(i,j);
                        }
                    }
                    if( state.needfgh )
                    {
                        for(j = 0; j <= n-1; j++)
                        {
                            for(k = 0; k <= n-1; k++)
                            {
                                state.h(j,k) = state.h(j,k)+2*a(i,j)*a(i,k);
                            }
                        }
                    }
                    if( state.needfij )
                    {
                        state.fi(i) = v-b(i);
                        ap::vmove(&state.j(i, 0), 1, &a(i, 0), 1, ap::vlen(0,n-1));
                    }
                }
                scerror = scerror||!rkindvsstatecheck(rkind, state);
            }
            minlmresults(state, x, rep);
            eqerror = eqerror||rep.terminationtype<=0;
            for(i = 0; i <= n-1; i++)
            {
                eqerror = eqerror||ap::fp_greater(fabs(x(i)-xe(i)),0.001);
            }
        }
    }
    
    //
    // Testing convergence properties using
    // different optimizer types and different conditions
    //
    s = 100;
    for(rkind = 0; rkind <= 2; rkind++)
    {
        for(ckind = 0; ckind <= 3; ckind++)
        {
            epsg = 0;
            epsf = 0;
            epsx = 0;
            maxits = 0;
            if( ckind==0 )
            {
                epsf = 0.0001;
            }
            if( ckind==1 )
            {
                epsx = 0.0001;
            }
            if( ckind==2 )
            {
                maxits = 2;
            }
            if( ckind==3 )
            {
                epsg = 0.0001;
            }
            x.setlength(3);
            n = 3;
            m = 3;
            for(i = 0; i <= 2; i++)
            {
                x(i) = 6;
            }
            if( rkind==0 )
            {
                minlmcreatefj(n, m, x, state);
            }
            if( rkind==1 )
            {
                minlmcreatefgj(n, m, x, state);
            }
            if( rkind==2 )
            {
                minlmcreatefgh(n, x, state);
            }
            minlmsetcond(state, epsg, epsf, epsx, maxits);
            while(minlmiteration(state))
            {
                if( state.needf||state.needfg||state.needfgh )
                {
                    state.f = s*ap::sqr(exp(state.x(0))-2)+ap::sqr(ap::sqr(state.x(1))+1)+ap::sqr(state.x(2)-state.x(0));
                }
                if( state.needfg||state.needfgh )
                {
                    state.g(0) = s*2*(exp(state.x(0))-2)*exp(state.x(0))+2*(state.x(0)-state.x(2));
                    state.g(1) = 2*(ap::sqr(state.x(1))+1)*2*state.x(1);
                    state.g(2) = 2*(state.x(2)-state.x(0));
                }
                if( state.needfgh )
                {
                    state.h(0,0) = s*(4*ap::sqr(exp(state.x(0)))-4*exp(state.x(0)))+2;
                    state.h(0,1) = 0;
                    state.h(0,2) = -2;
                    state.h(1,0) = 0;
                    state.h(1,1) = 12*ap::sqr(state.x(1))+4;
                    state.h(1,2) = 0;
                    state.h(2,0) = -2;
                    state.h(2,1) = 0;
                    state.h(2,2) = 2;
                }
                if( state.needfij )
                {
                    state.fi(0) = s*(exp(state.x(0))-2);
                    state.j(0,0) = s*exp(state.x(0));
                    state.j(0,1) = 0;
                    state.j(0,2) = 0;
                    state.fi(1) = ap::sqr(state.x(1))+1;
                    state.j(1,0) = 0;
                    state.j(1,1) = 2*state.x(1);
                    state.j(1,2) = 0;
                    state.fi(2) = state.x(2)-state.x(0);
                    state.j(2,0) = -1;
                    state.j(2,1) = 0;
                    state.j(2,2) = 1;
                }
                scerror = scerror||!rkindvsstatecheck(rkind, state);
            }
            minlmresults(state, x, rep);
            if( ckind==0 )
            {
                converror = converror||ap::fp_greater(fabs(x(0)-log(double(2))),0.05);
                converror = converror||ap::fp_greater(fabs(x(1)),0.05);
                converror = converror||ap::fp_greater(fabs(x(2)-log(double(2))),0.05);
                converror = converror||rep.terminationtype!=1;
            }
            if( ckind==1 )
            {
                converror = converror||ap::fp_greater(fabs(x(0)-log(double(2))),0.05);
                converror = converror||ap::fp_greater(fabs(x(1)),0.05);
                converror = converror||ap::fp_greater(fabs(x(2)-log(double(2))),0.05);
                converror = converror||rep.terminationtype!=2;
            }
            if( ckind==2 )
            {
                converror = converror||rep.terminationtype!=5||rep.iterationscount!=maxits;
            }
            if( ckind==3 )
            {
                converror = converror||ap::fp_greater(fabs(x(0)-log(double(2))),0.05);
                converror = converror||ap::fp_greater(fabs(x(1)),0.05);
                converror = converror||ap::fp_greater(fabs(x(2)-log(double(2))),0.05);
                converror = converror||rep.terminationtype!=4;
            }
        }
    }
    
    //
    // Other properties:
    // 1. test reports (F should form monotone sequence)
    // 2. test maximum step
    //
    for(rkind = 0; rkind <= 2; rkind++)
    {
        
        //
        // reports:
        // * check that first report is initial point
        // * check that F is monotone decreasing
        // * check that last report is final result
        //
        n = 3;
        m = 3;
        s = 100;
        x.setlength(n);
        xlast.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            x(i) = 6;
        }
        if( rkind==0 )
        {
            minlmcreatefj(n, m, x, state);
        }
        if( rkind==1 )
        {
            minlmcreatefgj(n, m, x, state);
        }
        if( rkind==2 )
        {
            minlmcreatefgh(n, x, state);
        }
        minlmsetcond(state, double(0), double(0), double(0), 4);
        minlmsetxrep(state, true);
        fprev = ap::maxrealnumber;
        while(minlmiteration(state))
        {
            if( state.needf||state.needfg||state.needfgh )
            {
                state.f = s*ap::sqr(exp(state.x(0))-2)+ap::sqr(state.x(1))+ap::sqr(state.x(2)-state.x(0));
            }
            if( state.needfg||state.needfgh )
            {
                state.g(0) = s*2*(exp(state.x(0))-2)*exp(state.x(0))+2*(state.x(0)-state.x(2));
                state.g(1) = 2*state.x(1);
                state.g(2) = 2*(state.x(2)-state.x(0));
            }
            if( state.needfgh )
            {
                state.h(0,0) = s*(4*ap::sqr(exp(state.x(0)))-4*exp(state.x(0)))+2;
                state.h(0,1) = 0;
                state.h(0,2) = -2;
                state.h(1,0) = 0;
                state.h(1,1) = 2;
                state.h(1,2) = 0;
                state.h(2,0) = -2;
                state.h(2,1) = 0;
                state.h(2,2) = 2;
            }
            if( state.needfij )
            {
                state.fi(0) = sqrt(s)*(exp(state.x(0))-2);
                state.j(0,0) = sqrt(s)*exp(state.x(0));
                state.j(0,1) = 0;
                state.j(0,2) = 0;
                state.fi(1) = state.x(1);
                state.j(1,0) = 0;
                state.j(1,1) = 1;
                state.j(1,2) = 0;
                state.fi(2) = state.x(2)-state.x(0);
                state.j(2,0) = -1;
                state.j(2,1) = 0;
                state.j(2,2) = 1;
            }
            scerror = scerror||!rkindvsstatecheck(rkind, state);
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
        minlmresults(state, x, rep);
        for(i = 0; i <= n-1; i++)
        {
            othererrors = othererrors||ap::fp_neq(x(i),xlast(i));
        }
    }
    n = 1;
    x.setlength(n);
    x(0) = 100;
    stpmax = 0.05+0.05*ap::randomreal();
    minlmcreatefgh(n, x, state);
    minlmsetcond(state, 1.0E-9, double(0), double(0), 0);
    minlmsetstpmax(state, stpmax);
    minlmsetxrep(state, true);
    xprev = x(0);
    while(minlmiteration(state))
    {
        if( state.needf||state.needfg||state.needfgh )
        {
            state.f = exp(state.x(0))+exp(-state.x(0));
        }
        if( state.needfg||state.needfgh )
        {
            state.g(0) = exp(state.x(0))-exp(-state.x(0));
        }
        if( state.needfgh )
        {
            state.h(0,0) = exp(state.x(0))+exp(-state.x(0));
        }
        othererrors = othererrors||ap::fp_greater(fabs(state.x(0)-xprev),(1+sqrt(ap::machineepsilon))*stpmax);
        if( state.xupdated )
        {
            xprev = state.x(0);
        }
    }
    
    //
    // end
    //
    waserrors = referror||lin1error||lin2error||eqerror||converror||scerror||othererrors;
    if( !silent )
    {
        printf("TESTING LEVENBERG-MARQUARDT OPTIMIZATION\n");
        printf("REFERENCE PROBLEM:                        ");
        if( referror )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("1-D PROBLEM #1:                           ");
        if( lin1error )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("1-D PROBLEM #2:                           ");
        if( lin2error )
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
        printf("STATE FIELDS CONSISTENCY:                 ");
        if( scerror )
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
Asserts that State fields are consistent with RKind.
Returns False otherwise.
*************************************************************************/
static bool rkindvsstatecheck(int rkind, const minlmstate& state)
{
    bool result;
    int nset;

    nset = 0;
    if( state.needf )
    {
        nset = nset+1;
    }
    if( state.needfg )
    {
        nset = nset+1;
    }
    if( state.needfij )
    {
        nset = nset+1;
    }
    if( state.needfgh )
    {
        nset = nset+1;
    }
    if( state.xupdated )
    {
        nset = nset+1;
    }
    if( nset!=1 )
    {
        result = false;
        return result;
    }
    if( rkind==0&&(state.needfg||state.needfgh) )
    {
        result = false;
        return result;
    }
    if( rkind==1&&state.needfgh )
    {
        result = false;
        return result;
    }
    if( rkind==2&&state.needfij )
    {
        result = false;
        return result;
    }
    result = true;
    return result;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testminlmunit_test_silent()
{
    bool result;

    result = testminlm(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testminlmunit_test()
{
    bool result;

    result = testminlm(false);
    return result;
}





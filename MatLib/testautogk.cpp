
 
#include <stdio.h>
#include "testautogk.h"

/*************************************************************************
Test
*************************************************************************/
bool testautogkunit(bool silent)
{
    bool result;
    double a;
    double b;
    autogkstate state;
    autogkreport rep;
    double v;
    double exact;
    double eabs;
    double alpha;
    int pkind;
    double errtol;
    bool simpleerrors;
    bool sngenderrors;
    bool waserrors;

    simpleerrors = false;
    sngenderrors = false;
    waserrors = false;
    errtol = 10000*ap::machineepsilon;
    
    //
    // Simple test: integral(exp(x),+-1,+-2), no maximum width requirements
    //
    a = (2*ap::randominteger(2)-1)*1.0;
    b = (2*ap::randominteger(2)-1)*2.0;
    autogksmooth(a, b, state);
    while(autogkiteration(state))
    {
        state.f = exp(state.x);
    }
    autogkresults(state, v, rep);
    exact = exp(b)-exp(a);
    eabs = fabs(exp(b)-exp(a));
    if( rep.terminationtype<=0 )
    {
        simpleerrors = true;
    }
    else
    {
        simpleerrors = simpleerrors||ap::fp_greater(fabs(exact-v),errtol*eabs);
    }
    
    //
    // Simple test: integral(exp(x),+-1,+-2), XWidth=0.1
    //
    a = (2*ap::randominteger(2)-1)*1.0;
    b = (2*ap::randominteger(2)-1)*2.0;
    autogksmoothw(a, b, 0.1, state);
    while(autogkiteration(state))
    {
        state.f = exp(state.x);
    }
    autogkresults(state, v, rep);
    exact = exp(b)-exp(a);
    eabs = fabs(exp(b)-exp(a));
    if( rep.terminationtype<=0 )
    {
        simpleerrors = true;
    }
    else
    {
        simpleerrors = simpleerrors||ap::fp_greater(fabs(exact-v),errtol*eabs);
    }
    
    //
    // Simple test: integral(cos(100*x),0,2*pi), no maximum width requirements
    //
    a = 0;
    b = 2*ap::pi();
    autogksmooth(a, b, state);
    while(autogkiteration(state))
    {
        state.f = cos(100*state.x);
    }
    autogkresults(state, v, rep);
    exact = 0;
    eabs = 4;
    if( rep.terminationtype<=0 )
    {
        simpleerrors = true;
    }
    else
    {
        simpleerrors = simpleerrors||ap::fp_greater(fabs(exact-v),errtol*eabs);
    }
    
    //
    // Simple test: integral(cos(100*x),0,2*pi), XWidth=0.3
    //
    a = 0;
    b = 2*ap::pi();
    autogksmoothw(a, b, 0.3, state);
    while(autogkiteration(state))
    {
        state.f = cos(100*state.x);
    }
    autogkresults(state, v, rep);
    exact = 0;
    eabs = 4;
    if( rep.terminationtype<=0 )
    {
        simpleerrors = true;
    }
    else
    {
        simpleerrors = simpleerrors||ap::fp_greater(fabs(exact-v),errtol*eabs);
    }
    
    //
    // singular problem on [a,b] = [0.1, 0.5]
    //     f2(x) = (1+x)*(b-x)^alpha, -1 < alpha < 1
    //
    for(pkind = 0; pkind <= 6; pkind++)
    {
        a = 0.1;
        b = 0.5;
        if( pkind==0 )
        {
            alpha = -0.9;
        }
        if( pkind==1 )
        {
            alpha = -0.5;
        }
        if( pkind==2 )
        {
            alpha = -0.1;
        }
        if( pkind==3 )
        {
            alpha = 0.0;
        }
        if( pkind==4 )
        {
            alpha = 0.1;
        }
        if( pkind==5 )
        {
            alpha = 0.5;
        }
        if( pkind==6 )
        {
            alpha = 0.9;
        }
        
        //
        // f1(x) = (1+x)*(x-a)^alpha, -1 < alpha < 1
        // 1. use singular integrator for [a,b]
        // 2. use singular integrator for [b,a]
        //
        exact = pow(b-a, alpha+2)/(alpha+2)+(1+a)*pow(b-a, alpha+1)/(alpha+1);
        eabs = fabs(exact);
        autogksingular(a, b, alpha, 0.0, state);
        while(autogkiteration(state))
        {
            if( ap::fp_less(state.xminusa,0.01) )
            {
                state.f = pow(state.xminusa, alpha)*(1+state.x);
            }
            else
            {
                state.f = pow(state.x-a, alpha)*(1+state.x);
            }
        }
        autogkresults(state, v, rep);
        if( rep.terminationtype<=0 )
        {
            sngenderrors = true;
        }
        else
        {
            sngenderrors = sngenderrors||ap::fp_greater(fabs(v-exact),errtol*eabs);
        }
        autogksingular(b, a, 0.0, alpha, state);
        while(autogkiteration(state))
        {
            if( ap::fp_greater(state.bminusx,-0.01) )
            {
                state.f = pow(-state.bminusx, alpha)*(1+state.x);
            }
            else
            {
                state.f = pow(state.x-a, alpha)*(1+state.x);
            }
        }
        autogkresults(state, v, rep);
        if( rep.terminationtype<=0 )
        {
            sngenderrors = true;
        }
        else
        {
            sngenderrors = sngenderrors||ap::fp_greater(fabs(-v-exact),errtol*eabs);
        }
        
        //
        // f1(x) = (1+x)*(b-x)^alpha, -1 < alpha < 1
        // 1. use singular integrator for [a,b]
        // 2. use singular integrator for [b,a]
        //
        exact = (1+b)*pow(b-a, alpha+1)/(alpha+1)-pow(b-a, alpha+2)/(alpha+2);
        eabs = fabs(exact);
        autogksingular(a, b, 0.0, alpha, state);
        while(autogkiteration(state))
        {
            if( ap::fp_less(state.bminusx,0.01) )
            {
                state.f = pow(state.bminusx, alpha)*(1+state.x);
            }
            else
            {
                state.f = pow(b-state.x, alpha)*(1+state.x);
            }
        }
        autogkresults(state, v, rep);
        if( rep.terminationtype<=0 )
        {
            sngenderrors = true;
        }
        else
        {
            sngenderrors = sngenderrors||ap::fp_greater(fabs(v-exact),errtol*eabs);
        }
        autogksingular(b, a, alpha, 0.0, state);
        while(autogkiteration(state))
        {
            if( ap::fp_greater(state.xminusa,-0.01) )
            {
                state.f = pow(-state.xminusa, alpha)*(1+state.x);
            }
            else
            {
                state.f = pow(b-state.x, alpha)*(1+state.x);
            }
        }
        autogkresults(state, v, rep);
        if( rep.terminationtype<=0 )
        {
            sngenderrors = true;
        }
        else
        {
            sngenderrors = sngenderrors||ap::fp_greater(fabs(-v-exact),errtol*eabs);
        }
    }
    
    //
    // end
    //
    waserrors = simpleerrors||sngenderrors;
    if( !silent )
    {
        printf("TESTING AUTOGK\n");
        printf("INTEGRATION WITH GIVEN ACCURACY:          ");
        if( simpleerrors||sngenderrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* SIMPLE PROBLEMS:                        ");
        if( simpleerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* SINGULAR PROBLEMS (ENDS OF INTERVAL):   ");
        if( sngenderrors )
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
Silent unit test
*************************************************************************/
bool testautogk_test_silent()
{
    bool result;

    result = testautogkunit(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testautogk_test()
{
    bool result;

    result = testautogkunit(false);
    return result;
}





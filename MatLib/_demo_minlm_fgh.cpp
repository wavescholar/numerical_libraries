#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "minlm.h"

int main(int argc, char **argv)
{
    minlmstate state;
    minlmreport rep;
    ap::real_1d_array s;
    double x;
    double y;

    
    //
    // Example of solving simple task using FGH scheme.
    //
    // Function minimized:
    //     F = (x-2*y)^2 + (x-2)^2 + (y-1)^2
    // exact solution is (2,1).
    //
    s.setlength(2);
    s(0) = ap::randomreal()-0.5;
    s(1) = ap::randomreal()-0.5;
    minlmcreatefgh(2, s, state);
    minlmsetcond(state, 0.0, 0.0, 0.001, 0);
    while(minlmiteration(state))
    {
        x = state.x(0);
        y = state.x(1);
        if( state.needf )
        {
            state.f = ap::sqr(x-2*y)+ap::sqr(x-2)+ap::sqr(y-1);
        }
        if( state.needfg )
        {
            state.f = ap::sqr(x-2*y)+ap::sqr(x-2)+ap::sqr(y-1);
            state.g(0) = 2*(x-2*y)+2*(x-2)+0;
            state.g(1) = -4*(x-2*y)+0+2*(y-1);
        }
        if( state.needfgh )
        {
            state.f = ap::sqr(x-2*y)+ap::sqr(x-2)+ap::sqr(y-1);
            state.g(0) = 2*(x-2*y)+2*(x-2)+0;
            state.g(1) = -4*(x-2*y)+0+2*(y-1);
            state.h(0,0) = 4;
            state.h(1,0) = -4;
            state.h(0,1) = -4;
            state.h(1,1) = 10;
        }
    }
    minlmresults(state, s, rep);
    
    //
    // output results
    //
    printf("X = %4.2lf (correct value - 2.00)\n",
        double(s(0)));
    printf("Y = %4.2lf (correct value - 1.00)\n",
        double(s(1)));
    printf("TerminationType = %0ld (should be 2 - stopping when step is small enough)\n",
        long(rep.terminationtype));
    printf("NFunc = %0ld\n",
        long(rep.nfunc));
    printf("NJac  = %0ld\n",
        long(rep.njac));
    printf("NGrad = %0ld\n",
        long(rep.ngrad));
    printf("NHess = %0ld\n",
        long(rep.nhess));
    return 0;
}


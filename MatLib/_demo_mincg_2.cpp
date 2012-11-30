#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "mincg.h"

int main(int argc, char **argv)
{
    int n;
    mincgstate state;
    mincgreport rep;
    ap::real_1d_array s;
    double x;
    double y;

    
    //
    // Function minimized:
    //     F = exp(x-1) + exp(1-x) + (y-x)^2
    // N = 2 - task dimension.
    //
    // Take a look at MinCGSetStpMax() call - it prevents us
    // from overflow (which may be result of too large step).
    // Try to comment it and see what will happen.
    //
    n = 2;
    s.setlength(2);
    s(0) = 10;
    s(1) = ap::randomreal()-0.5;
    mincgcreate(n, s, state);
    mincgsetcond(state, 0.0, 0.0, 0.0001, 0);
    mincgsetxrep(state, true);
    mincgsetstpmax(state, 1.0);
    printf("\n\nF = exp(x-1) + exp(1-x) + (y-x)^2\n");
    printf("OPTIMIZATION STARTED\n");
    while(mincgiteration(state))
    {
        if( state.needfg )
        {
            x = state.x(0);
            y = state.x(1);
            state.f = exp(x-1)+exp(1-x)+ap::sqr(y-x);
            state.g(0) = exp(x-1)-exp(1-x)+2*(x-y);
            state.g(1) = 2*(y-x);
        }
        if( state.xupdated )
        {
            printf("    F(%8.5lf,%8.5lf)=%0.5lf\n",
                double(state.x(0)),
                double(state.x(1)),
                double(state.f));
        }
    }
    printf("OPTIMIZATION STOPPED\n");
    mincgresults(state, s, rep);
    
    //
    // output results
    //
    printf("X = %4.2lf (should be 1.00)\n",
        double(s(0)));
    printf("Y = %4.2lf (should be 1.00)\n\n\n",
        double(s(1)));
    return 0;
}


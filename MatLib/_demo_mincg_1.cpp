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
    //     F = (x-1)^4 + (y-x)^2
    // N = 2 - task dimension.
    //
    n = 2;
    s.setlength(2);
    s(0) = 10;
    s(1) = 11;
    mincgcreate(n, s, state);
    mincgsetcond(state, 0.0, 0.0, 0.00001, 0);
    mincgsetxrep(state, true);
    printf("\n\nF = (x-1)^4 + (y-x)^2\n");
    printf("OPTIMIZATION STARTED\n");
    while(mincgiteration(state))
    {
        if( state.needfg )
        {
            x = state.x(0);
            y = state.x(1);
            state.f = ap::sqr(ap::sqr(x-1))+ap::sqr(y-x);
            state.g(0) = 4*ap::sqr(x-1)*(x-1)+2*(x-y);
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


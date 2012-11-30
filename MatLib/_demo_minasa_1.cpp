#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "minasa.h"

int main(int argc, char **argv)
{
    int n;
    int i;
    minasastate state;
    minasareport rep;
    ap::real_1d_array s;
    ap::real_1d_array bndl;
    ap::real_1d_array bndu;
    double x;
    double y;
    double z;

    
    //
    // Function being minimized:
    //     F = x+2y+3z subject to 0<=x<=1, 0<=y<=1, 0<=z<=1.
    //
    n = 3;
    s.setlength(n);
    bndl.setlength(n);
    bndu.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        s(i) = 1;
        bndl(i) = 0;
        bndu(i) = 1;
    }
    minasacreate(n, s, bndl, bndu, state);
    minasasetcond(state, 0.0, 0.0, 0.00001, 0);
    minasasetxrep(state, true);
    printf("\n\nF = x+2y+3z subject to 0<=x<=1, 0<=y<=1, 0<=z<=1\n");
    printf("OPTIMIZATION STARTED\n");
    while(minasaiteration(state))
    {
        if( state.needfg )
        {
            x = state.x(0);
            y = state.x(1);
            z = state.x(2);
            state.f = x+2*y+3*z;
            state.g(0) = 1;
            state.g(1) = 2;
            state.g(2) = 3;
        }
        if( state.xupdated )
        {
            printf("    F(%4.2lf,%4.2lf,%4.2lf)=%0.3lf\n",
                double(state.x(0)),
                double(state.x(1)),
                double(state.x(2)),
                double(state.f));
        }
    }
    printf("OPTIMIZATION STOPPED\n");
    minasaresults(state, s, rep);
    
    //
    // output results
    //
    printf("X = %4.2lf (should be 0.00)\n",
        double(s(0)));
    printf("Y = %4.2lf (should be 0.00)\n",
        double(s(1)));
    printf("Z = %4.2lf (should be 0.00)\n\n\n",
        double(s(2)));
    return 0;
}


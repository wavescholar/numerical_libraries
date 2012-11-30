#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "odesolver.h"

int main(int argc, char **argv)
{
    ap::real_1d_array x;
    ap::real_1d_array y;
    ap::real_2d_array ytbl;
    double eps;
    double h;
    int m;
    int i;
    odesolverstate state;
    odesolverreport rep;

    
    //
    // ODESolver unit is used to solve simple ODE:
    // y' = y, y(0) = 1.
    //
    // Its solution is well known in academic circles :)
    //
    // No intermediate values are calculated,
    // just starting and final points.
    //
    y.setlength(1);
    y(0) = 1;
    x.setlength(2);
    x(0) = 0;
    x(1) = 1;
    eps = 1.0E-4;
    h = 0.01;
    odesolverrkck(y, 1, x, 2, eps, h, state);
    while(odesolveriteration(state))
    {
        state.dy(0) = state.y(0);
    }
    odesolverresults(state, m, x, ytbl, rep);
    printf("    X  Y(X)\n");
    for(i = 0; i <= m-1; i++)
    {
        printf("%5.3lf %5.3lf\n",
            double(x(i)),
            double(ytbl(i,0)));
    }
    return 0;
}


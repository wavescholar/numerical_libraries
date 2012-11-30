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
    // y'' = -y, y(0) = 0, y'(0)=1.
    //
    // This ODE may be written as first-order system:
    // y' =  z
    // z' = -y
    //
    // Its solution is well known in academic circles :)
    //
    // Three intermediate values are calculated,
    // plus starting and final points.
    //
    y.setlength(2);
    y(0) = 0;
    y(1) = 1;
    x.setlength(5);
    x(0) = ap::pi()*0/4;
    x(1) = ap::pi()*1/4;
    x(2) = ap::pi()*2/4;
    x(3) = ap::pi()*3/4;
    x(4) = ap::pi()*4/4;
    eps = 1.0E-8;
    h = 0.01;
    odesolverrkck(y, 2, x, 5, eps, h, state);
    while(odesolveriteration(state))
    {
        state.dy(0) = state.y(1);
        state.dy(1) = -state.y(0);
    }
    odesolverresults(state, m, x, ytbl, rep);
    printf("     X   Y(X)     Error\n");
    for(i = 0; i <= m-1; i++)
    {
        printf("%6.3lf %6.3lf  %8.1le\n",
            double(x(i)),
            double(ytbl(i,0)),
            double(fabs(ytbl(i,0)-sin(x(i)))));
    }
    return 0;
}


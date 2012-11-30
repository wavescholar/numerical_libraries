#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "autogk.h"

int main(int argc, char **argv)
{
    autogkstate state;
    double v;
    autogkreport rep;
    double a;
    double b;
    double alpha;

    
    //
    // f1(x) = (1+x)*(x-a)^alpha, alpha=-0.3
    // Exact answer is (B-A)^(Alpha+2)/(Alpha+2) + (1+A)*(B-A)^(Alpha+1)/(Alpha+1)
    //
    // This code demonstrates use of the State.XMinusA (State.BMinusX) field.
    //
    // If we try to use State.X instead of State.XMinusA,
    // we will end up dividing by zero! (in 64-bit precision)
    //
    a = 1.0;
    b = 5.0;
    alpha = -0.9;
    autogksingular(a, b, alpha, 0.0, state);
    while(autogkiteration(state))
    {
        state.f = pow(state.xminusa, alpha)*(1+state.x);
    }
    autogkresults(state, v, rep);
    printf("integral((1+x)*(x-a)^alpha) on [%0.1lf; %0.1lf] = %0.2lf\nExact answer is %0.2lf\n",
        double(a),
        double(b),
        double(v),
        double(pow(b-a, alpha+2)/(alpha+2)+(1+a)*pow(b-a, alpha+1)/(alpha+1)));
    return 0;
}


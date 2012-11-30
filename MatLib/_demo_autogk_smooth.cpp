#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "autogk.h"

int main(int argc, char **argv)
{
    autogkstate state;
    double v;
    autogkreport rep;

    
    //
    // f(x) = x*sin(x), integrated at [-pi, pi].
    // Exact answer is 2*pi
    //
    autogksmooth(-ap::pi(), +ap::pi(), state);
    while(autogkiteration(state))
    {
        state.f = state.x*sin(state.x);
    }
    autogkresults(state, v, rep);
    printf("integral(x*sin(x),-pi,+pi) = %0.2lf\nExact answer is %0.2lf\n",
        double(v),
        double(2*ap::pi()));
    return 0;
}


#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "minlm.h"

int main(int argc, char **argv)
{
    minlmstate state;
    minlmreport rep;
    int i;
    ap::real_1d_array s;
    ap::real_1d_array x;
    ap::real_1d_array y;
    double fi;
    int n;
    int m;

    
    //
    // Example of solving polynomial approximation task using FJ scheme.
    //
    // Data points:
    //     xi are random numbers from [-1,+1],
    //
    // Function being fitted:
    //     yi = exp(xi) - sin(xi) - x^3/3
    //
    // Function being minimized:
    //     F(a,b,c) =
    //         (a + b*x0 + c*x0^2 - y0)^2 +
    //         (a + b*x1 + c*x1^2 - y1)^2 + ...
    //
    n = 3;
    s.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        s(i) = ap::randomreal()-0.5;
    }
    m = 100;
    x.setlength(m);
    y.setlength(m);
    for(i = 0; i <= m-1; i++)
    {
        x(i) = double(2*i)/double(m-1)-1;
        y(i) = exp(x(i))-sin(x(i))-x(i)*x(i)*x(i)/3;
    }
    
    //
    // Now S stores starting point, X and Y store points being fitted.
    //
    minlmcreatefj(n, m, s, state);
    minlmsetcond(state, 0.0, 0.0, 0.001, 0);
    while(minlmiteration(state))
    {
        if( state.needf )
        {
            state.f = 0;
        }
        for(i = 0; i <= m-1; i++)
        {
            
            //
            // "a" is stored in State.X[0]
            // "b" - State.X[1]
            // "c" - State.X[2]
            //
            fi = state.x(0)+state.x(1)*x(i)+state.x(2)*ap::sqr(x(i))-y(i);
            if( state.needf )
            {
                
                //
                // F is equal to sum of fi squared.
                //
                state.f = state.f+ap::sqr(fi);
            }
            if( state.needfij )
            {
                
                //
                // Fi
                //
                state.fi(i) = fi;
                
                //
                // dFi/da
                //
                state.j(i,0) = 1;
                
                //
                // dFi/db
                //
                state.j(i,1) = x(i);
                
                //
                // dFi/dc
                //
                state.j(i,2) = ap::sqr(x(i));
            }
        }
    }
    minlmresults(state, s, rep);
    
    //
    // output results
    //
    printf("A = %4.2lf\n",
        double(s(0)));
    printf("B = %4.2lf\n",
        double(s(1)));
    printf("C = %4.2lf\n",
        double(s(2)));
    printf("TerminationType = %0ld (should be 2 - stopping when step is small enough)\n",
        long(rep.terminationtype));
    return 0;
}


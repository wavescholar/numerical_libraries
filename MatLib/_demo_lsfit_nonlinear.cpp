#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "lsfit.h"

int main(int argc, char **argv)
{
    int m;
    int n;
    int k;
    ap::real_1d_array y;
    ap::real_2d_array x;
    ap::real_1d_array c;
    lsfitreport rep;
    lsfitstate state;
    int info;
    double epsf;
    double epsx;
    int maxits;
    int i;
    int j;
    double a;
    double b;

    printf("Fitting 0.5(1+cos(x)) on [-pi,+pi] with exp(-alpha*x^2)\n");
    
    //
    // Fitting 0.5(1+cos(x)) on [-pi,+pi] with Gaussian exp(-alpha*x^2):
    // * without Hessian (gradient only)
    // * using alpha=1 as initial value
    // * using 1000 uniformly distributed points to fit to
    //
    // Notes:
    // * N - number of points
    // * M - dimension of space where points reside
    // * K - number of parameters being fitted
    //
    n = 1000;
    m = 1;
    k = 1;
    a = -ap::pi();
    b = +ap::pi();
    
    //
    // Prepare task matrix
    //
    y.setlength(n);
    x.setlength(n, m);
    c.setlength(k);
    for(i = 0; i <= n-1; i++)
    {
        x(i,0) = a+(b-a)*i/(n-1);
        y(i) = 0.5*(1+cos(x(i,0)));
    }
    c(0) = 1.0;
    epsf = 0.0;
    epsx = 0.0001;
    maxits = 0;
    
    //
    // Solve
    //
    lsfitnonlinearfg(x, y, c, n, m, k, true, state);
    lsfitnonlinearsetcond(state, epsf, epsx, maxits);
    while(lsfitnonlineariteration(state))
    {
        if( state.needf )
        {
            
            //
            // F(x) = Exp(-alpha*x^2)
            //
            state.f = exp(-state.c(0)*ap::sqr(state.x(0)));
        }
        if( state.needfg )
        {
            
            //
            // F(x)      = Exp(-alpha*x^2)
            // dF/dAlpha = (-x^2)*Exp(-alpha*x^2)
            //
            state.f = exp(-state.c(0)*ap::sqr(state.x(0)));
            state.g(0) = -ap::sqr(state.x(0))*state.f;
        }
    }
    lsfitnonlinearresults(state, info, c, rep);
    printf("alpha:   %0.3lf\n",
        double(c(0)));
    printf("rms.err: %0.3lf\n",
        double(rep.rmserror));
    printf("max.err: %0.3lf\n",
        double(rep.maxerror));
    printf("Termination type: %0ld\n",
        long(info));
    printf("\n\n");
    return 0;
}


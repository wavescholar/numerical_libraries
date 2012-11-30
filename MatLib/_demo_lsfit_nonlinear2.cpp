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

    printf("Fitting 1-x^2 on [-1,+1] with cos(alpha*pi*x)+beta\n");
    
    //
    // Fitting 1-x^2 on [-1,+1] with cos(alpha*pi*x)+beta:
    // * using Hessian
    // * using alpha=1 and beta=0 as initial values
    // * using 1000 uniformly distributed points to fit to
    //
    // Notes:
    // * N - number of points
    // * M - dimension of space where points reside
    // * K - number of parameters being fitted
    //
    n = 1000;
    m = 1;
    k = 2;
    a = -1;
    b = +1;
    
    //
    // Prepare task matrix
    //
    y.setlength(n);
    x.setlength(n, m);
    c.setlength(k);
    for(i = 0; i <= n-1; i++)
    {
        x(i,0) = a+(b-a)*i/(n-1);
        y(i) = 1-ap::sqr(x(i,0));
    }
    c(0) = 1.0;
    c(1) = 0.0;
    epsf = 0.0;
    epsx = 0.0001;
    maxits = 0;
    
    //
    // Solve
    //
    lsfitnonlinearfgh(x, y, c, n, m, k, state);
    lsfitnonlinearsetcond(state, epsf, epsx, maxits);
    while(lsfitnonlineariteration(state))
    {
        
        //
        // F(x) = Cos(alpha*pi*x)+beta
        //
        state.f = cos(state.c(0)*ap::pi()*state.x(0))+state.c(1);
        
        //
        // F(x)      = Cos(alpha*pi*x)+beta
        // dF/dAlpha = -pi*x*Sin(alpha*pi*x)
        // dF/dBeta  = 1.0
        //
        if( state.needfg||state.needfgh )
        {
            state.g(0) = -ap::pi()*state.x(0)*sin(state.c(0)*ap::pi()*state.x(0));
            state.g(1) = 1.0;
        }
        
        //
        // F(x)            = Cos(alpha*pi*x)+beta
        // d2F/dAlpha2     = -(pi*x)^2*Cos(alpha*pi*x)
        // d2F/dAlphadBeta = 0
        // d2F/dBeta2     =  0
        //
        if( state.needfgh )
        {
            state.h(0,0) = -ap::sqr(ap::pi()*state.x(0))*cos(state.c(0)*ap::pi()*state.x(0));
            state.h(0,1) = 0.0;
            state.h(1,0) = 0.0;
            state.h(1,1) = 0.0;
        }
    }
    lsfitnonlinearresults(state, info, c, rep);
    printf("alpha:   %0.3lf\n",
        double(c(0)));
    printf("beta:    %0.3lf\n",
        double(c(1)));
    printf("rms.err: %0.3lf\n",
        double(rep.rmserror));
    printf("max.err: %0.3lf\n",
        double(rep.maxerror));
    printf("Termination type: %0ld\n",
        long(info));
    printf("\n\n");
    return 0;
}


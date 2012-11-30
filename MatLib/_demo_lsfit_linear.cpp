#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "lsfit.h"

int main(int argc, char **argv)
{
    int m;
    int n;
    ap::real_1d_array y;
    ap::real_2d_array fmatrix;
    ap::real_2d_array cmatrix;
    lsfitreport rep;
    int info;
    ap::real_1d_array c;
    int i;
    int j;
    double x;
    double a;
    double b;

    printf("\n\nFitting tan(x) by third degree polynomial\n\n");
    printf("Fit type             rms.err max.err    p(0)   dp(0)\n");
    
    //
    // Fitting tan(x) at [0, 0.4*pi] by third degree polynomial:
    // a) without constraints
    // b) constrained at x=0: p(0)=0
    // c) constrained at x=0: p'(0)=1
    // c) constrained at x=0: p(0)=0, p'(0)=1
    //
    m = 4;
    n = 100;
    a = 0;
    b = 0.4*ap::pi();
    
    //
    // Prepare task matrix
    //
    y.setlength(n);
    fmatrix.setlength(n, m);
    for(i = 0; i <= n-1; i++)
    {
        x = a+(b-a)*i/(n-1);
        y(i) = tan(x);
        fmatrix(i,0) = 1.0;
        for(j = 1; j <= m-1; j++)
        {
            fmatrix(i,j) = x*fmatrix(i,j-1);
        }
    }
    
    //
    // Solve unconstrained task
    //
    lsfitlinear(y, fmatrix, n, m, info, c, rep);
    printf("Unconstrained        %7.4lf %7.4lf %7.4lf %7.4lf\n",
        double(rep.rmserror),
        double(rep.maxerror),
        double(c(0)),
        double(c(1)));
    
    //
    // Solve constrained task, p(0)=0
    // Prepare constraints matrix:
    // * first M columns store values of basis functions at X=0
    // * last column stores zero (desired value at X=0)
    //
    cmatrix.setlength(1, m+1);
    cmatrix(0,0) = 1;
    for(i = 1; i <= m-1; i++)
    {
        cmatrix(0,i) = 0;
    }
    cmatrix(0,m) = 0;
    lsfitlinearc(y, fmatrix, cmatrix, n, m, 1, info, c, rep);
    printf("Constrained, p(0)=0  %7.4lf %7.4lf %7.4lf %7.4lf\n",
        double(rep.rmserror),
        double(rep.maxerror),
        double(c(0)),
        double(c(1)));
    
    //
    // Solve constrained task, p'(0)=0
    // Prepare constraints matrix:
    // * first M columns store derivatives of basis functions at X=0
    // * last column stores 1.0 (desired derivative at X=0)
    //
    cmatrix.setlength(1, m+1);
    for(i = 0; i <= m-1; i++)
    {
        cmatrix(0,i) = 0;
    }
    cmatrix(0,1) = 1;
    cmatrix(0,m) = 1;
    lsfitlinearc(y, fmatrix, cmatrix, n, m, 1, info, c, rep);
    printf("Constrained, dp(0)=1 %7.4lf %7.4lf %7.4lf %7.4lf\n",
        double(rep.rmserror),
        double(rep.maxerror),
        double(c(0)),
        double(c(1)));
    
    //
    // Solve constrained task, p(0)=0, p'(0)=0
    // Prepare constraints matrix:
    // * first M columns store values/derivatives of basis functions at X=0
    // * last column stores desired values/derivative at X=0
    //
    cmatrix.setlength(2, m+1);
    cmatrix(0,0) = 1;
    for(i = 1; i <= m-1; i++)
    {
        cmatrix(0,i) = 0;
    }
    cmatrix(0,m) = 0;
    for(i = 0; i <= m-1; i++)
    {
        cmatrix(1,i) = 0;
    }
    cmatrix(1,1) = 1;
    cmatrix(1,m) = 1;
    lsfitlinearc(y, fmatrix, cmatrix, n, m, 2, info, c, rep);
    printf("Constrained, both    %7.4lf %7.4lf %7.4lf %7.4lf\n",
        double(rep.rmserror),
        double(rep.maxerror),
        double(c(0)),
        double(c(1)));
    printf("\n\n");
    return 0;
}


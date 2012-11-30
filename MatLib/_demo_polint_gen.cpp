#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "polint.h"

int main(int argc, char **argv)
{
    ap::real_1d_array x;
    ap::real_1d_array y;
    int n;
    int i;
    double t;
    barycentricinterpolant p;
    double v;
    double dv;
    double d2v;
    double err;
    double maxerr;

    
    //
    // Demonstration
    //
    printf("POLYNOMIAL INTERPOLATION\n\n");
    printf("F(x)=sin(x), [0, pi]\n");
    printf("Second degree polynomial is used\n\n");
    
    //
    // Create polynomial interpolant
    //
    n = 3;
    x.setlength(n);
    y.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        x(i) = ap::pi()*i/(n-1);
        y(i) = sin(x(i));
    }
    polynomialbuild(x, y, n, p);
    
    //
    // Output results
    //
    barycentricdiff2(p, double(0), v, dv, d2v);
    printf("                 P(x)    F(x) \n");
    printf("function       %6.3lf  %6.3lf \n",
        double(barycentriccalc(p, double(0))),
        double(0));
    printf("d/dx(0)        %6.3lf  %6.3lf \n",
        double(dv),
        double(1));
    printf("d2/dx2(0)      %6.3lf  %6.3lf \n",
        double(d2v),
        double(0));
    printf("\n\n");
    return 0;
}


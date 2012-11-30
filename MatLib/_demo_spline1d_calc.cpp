#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "spline1d.h"

int main(int argc, char **argv)
{
    ap::real_1d_array x;
    ap::real_1d_array y;
    int n;
    int i;
    double t;
    spline1dinterpolant s;
    double v;
    double dv;
    double d2v;
    double err;
    double maxerr;

    
    //
    // Demonstration of Spline1DCalc(), Spline1DDiff(), Spline1DIntegrate()
    //
    printf("DEMONSTRATION OF Spline1DCalc(), Spline1DDiff(), Spline1DIntegrate()\n\n");
    printf("F(x)=sin(x), [0, pi]\n");
    printf("Natural cubic spline with 3 nodes is used\n\n");
    
    //
    // Create spline
    //
    n = 3;
    x.setlength(n);
    y.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        x(i) = ap::pi()*i/(n-1);
        y(i) = sin(x(i));
    }
    spline1dbuildcubic(x, y, n, 2, 0.0, 2, 0.0, s);
    
    //
    // Output results
    //
    spline1ddiff(s, double(0), v, dv, d2v);
    printf("                 S(x)    F(x) \n");
    printf("function       %6.3lf  %6.3lf \n",
        double(spline1dcalc(s, double(0))),
        double(0));
    printf("d/dx(0)        %6.3lf  %6.3lf \n",
        double(dv),
        double(1));
    printf("d2/dx2(0)      %6.3lf  %6.3lf \n",
        double(d2v),
        double(0));
    printf("integral(0,pi) %6.3lf  %6.3lf \n",
        double(spline1dintegrate(s, ap::pi())),
        double(2));
    printf("\n\n");
    return 0;
}


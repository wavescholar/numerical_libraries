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
    double err;
    double maxerr;

    
    //
    // Interpolation by natural Cubic spline.
    //
    printf("INTERPOLATION BY NATURAL CUBIC SPLINE\n\n");
    printf("F(x)=sin(x), [0, pi], 3 nodes\n\n");
    printf("     x   F(x)   S(x)  Error\n");
    
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
    spline1dbuildcubic(x, y, n, 1, double(+1), 1, double(-1), s);
    
    //
    // Output results
    //
    t = 0;
    maxerr = 0;
    while(ap::fp_less(t,0.999999*ap::pi()))
    {
        err = fabs(spline1dcalc(s, t)-sin(t));
        maxerr = ap::maxreal(err, maxerr);
        printf("%6.3lf %6.3lf %6.3lf %6.3lf\n",
            double(t),
            double(sin(t)),
            double(spline1dcalc(s, t)),
            double(err));
        t = ap::minreal(ap::pi(), t+0.25);
    }
    err = fabs(spline1dcalc(s, ap::pi())-sin(ap::pi()));
    maxerr = ap::maxreal(err, maxerr);
    printf("%6.3lf %6.3lf %6.3lf %6.3lf\n\n",
        double(ap::pi()),
        double(sin(ap::pi())),
        double(spline1dcalc(s, ap::pi())),
        double(err));
    printf("max|error| = %0.3lf\n",
        double(maxerr));
    printf("Try other demos (spline1d_linear, spline1d_hermite) and compare errors...\n\n\n");
    return 0;
}


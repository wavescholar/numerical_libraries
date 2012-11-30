#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "mlpbase.h"

int main(int argc, char **argv)
{
    multilayerperceptron net;
    ap::real_1d_array x;
    ap::real_1d_array y;

    
    //
    // regression task with 2 inputs (independent variables)
    // and 2 outputs (dependent variables).
    //
    // network weights are initialized with small random values.
    //
    mlpcreate0(2, 2, net);
    x.setlength(2);
    y.setlength(2);
    x(0) = ap::randomreal()-0.5;
    x(1) = ap::randomreal()-0.5;
    mlpprocess(net, x, y);
    printf("Regression task\n");
    printf("IN[0]  = %5.2lf\n",
        double(x(0)));
    printf("IN[1]  = %5.2lf\n",
        double(x(1)));
    printf("OUT[0] = %5.2lf\n",
        double(y(0)));
    printf("OUT[1] = %5.2lf\n",
        double(y(1)));
    return 0;
}


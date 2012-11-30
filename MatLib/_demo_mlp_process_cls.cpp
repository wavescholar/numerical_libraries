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
    // classification task with 2 inputs and 3 classes.
    //
    // network weights are initialized with small random values.
    //
    mlpcreatec0(2, 3, net);
    x.setlength(2);
    y.setlength(3);
    x(0) = ap::randomreal()-0.5;
    x(1) = ap::randomreal()-0.5;
    mlpprocess(net, x, y);
    
    //
    // output results
    //
    printf("Classification task\n");
    printf("IN[0]  = %5.2lf\n",
        double(x(0)));
    printf("IN[1]  = %5.2lf\n",
        double(x(1)));
    printf("Prob(Class=0|IN) = %5.2lf\n",
        double(y(0)));
    printf("Prob(Class=1|IN) = %5.2lf\n",
        double(y(1)));
    printf("Prob(Class=2|IN) = %5.2lf\n",
        double(y(2)));
    return 0;
}


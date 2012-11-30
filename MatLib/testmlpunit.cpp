
 
#include <stdio.h>
#include "testmlpunit.h"

static void createnetwork(multilayerperceptron& network,
     int nkind,
     double a1,
     double a2,
     int nin,
     int nhid1,
     int nhid2,
     int nout);
static void unsetnetwork(multilayerperceptron& network);
static void testinformational(int nkind,
     int nin,
     int nhid1,
     int nhid2,
     int nout,
     int passcount,
     bool& err);
static void testprocessing(int nkind,
     int nin,
     int nhid1,
     int nhid2,
     int nout,
     int passcount,
     bool& err);
static void testgradient(int nkind,
     int nin,
     int nhid1,
     int nhid2,
     int nout,
     int passcount,
     bool& err);
static void testhessian(int nkind,
     int nin,
     int nhid1,
     int nhid2,
     int nout,
     int passcount,
     bool& err);

bool testmlp(bool silent)
{
    bool result;
    bool waserrors;
    int passcount;
    int maxn;
    int maxhid;
    int info;
    int nf;
    int nhid;
    int nl;
    int nhid1;
    int nhid2;
    int nkind;
    int i;
    int j;
    multilayerperceptron network;
    multilayerperceptron network2;
    mlpreport rep;
    mlpcvreport cvrep;
    int ncount;
    ap::real_2d_array xy;
    ap::real_2d_array valxy;
    int ssize;
    int valsize;
    bool allsame;
    bool inferrors;
    bool procerrors;
    bool graderrors;
    bool hesserrors;
    bool trnerrors;

    waserrors = false;
    inferrors = false;
    procerrors = false;
    graderrors = false;
    hesserrors = false;
    trnerrors = false;
    passcount = 10;
    maxn = 4;
    maxhid = 4;
    
    //
    // General multilayer network tests
    //
    for(nf = 1; nf <= maxn; nf++)
    {
        for(nl = 1; nl <= maxn; nl++)
        {
            for(nhid1 = 0; nhid1 <= maxhid; nhid1++)
            {
                for(nhid2 = 0; nhid2 <= 0; nhid2++)
                {
                    for(nkind = 0; nkind <= 3; nkind++)
                    {
                        
                        //
                        //  Skip meaningless parameters combinations
                        //
                        if( nkind==1&&nl<2 )
                        {
                            continue;
                        }
                        if( nhid1==0&&nhid2!=0 )
                        {
                            continue;
                        }
                        
                        //
                        // Tests
                        //
                        testinformational(nkind, nf, nhid1, nhid2, nl, passcount, inferrors);
                        testprocessing(nkind, nf, nhid1, nhid2, nl, passcount, procerrors);
                        testgradient(nkind, nf, nhid1, nhid2, nl, passcount, graderrors);
                        testhessian(nkind, nf, nhid1, nhid2, nl, passcount, hesserrors);
                    }
                }
            }
        }
    }
    
    //
    // Test network training on simple XOR problem
    //
    xy.setbounds(0, 3, 0, 2);
    xy(0,0) = -1;
    xy(0,1) = -1;
    xy(0,2) = -1;
    xy(1,0) = +1;
    xy(1,1) = -1;
    xy(1,2) = +1;
    xy(2,0) = -1;
    xy(2,1) = +1;
    xy(2,2) = +1;
    xy(3,0) = +1;
    xy(3,1) = +1;
    xy(3,2) = -1;
    mlpcreate1(2, 2, 1, network);
    mlptrainlm(network, xy, 4, 0.001, 10, info, rep);
    trnerrors = trnerrors||ap::fp_greater(mlprmserror(network, xy, 4),0.1);
    
    //
    // Test CV on random noisy problem
    //
    ncount = 100;
    xy.setbounds(0, ncount-1, 0, 1);
    for(i = 0; i <= ncount-1; i++)
    {
        xy(i,0) = 2*ap::randomreal()-1;
        xy(i,1) = ap::randominteger(4);
    }
    mlpcreatec0(1, 4, network);
    mlpkfoldcvlm(network, xy, ncount, 0.001, 5, 10, info, rep, cvrep);
    
    //
    // Final report
    //
    waserrors = inferrors||procerrors||graderrors||hesserrors||trnerrors;
    if( !silent )
    {
        printf("MLP TEST\n");
        printf("INFORMATIONAL FUNCTIONS:                 ");
        if( !inferrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("BASIC PROCESSING:                        ");
        if( !procerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("GRADIENT CALCULATION:                    ");
        if( !graderrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("HESSIAN CALCULATION:                     ");
        if( !hesserrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("TRAINING:                                ");
        if( !trnerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        if( waserrors )
        {
            printf("TEST SUMMARY: FAILED\n");
        }
        else
        {
            printf("TEST SUMMARY: PASSED\n");
        }
        printf("\n\n");
    }
    result = !waserrors;
    return result;
}


/*************************************************************************
Network creation
*************************************************************************/
static void createnetwork(multilayerperceptron& network,
     int nkind,
     double a1,
     double a2,
     int nin,
     int nhid1,
     int nhid2,
     int nout)
{

    ap::ap_error::make_assertion(nin>0&&nhid1>=0&&nhid2>=0&&nout>0, "CreateNetwork error");
    ap::ap_error::make_assertion(nhid1!=0||nhid2==0, "CreateNetwork error");
    ap::ap_error::make_assertion(nkind!=1||nout>=2, "CreateNetwork error");
    if( nhid1==0 )
    {
        
        //
        // No hidden layers
        //
        if( nkind==0 )
        {
            mlpcreate0(nin, nout, network);
        }
        else
        {
            if( nkind==1 )
            {
                mlpcreatec0(nin, nout, network);
            }
            else
            {
                if( nkind==2 )
                {
                    mlpcreateb0(nin, nout, a1, a2, network);
                }
                else
                {
                    if( nkind==3 )
                    {
                        mlpcreater0(nin, nout, a1, a2, network);
                    }
                }
            }
        }
        return;
    }
    if( nhid2==0 )
    {
        
        //
        // One hidden layer
        //
        if( nkind==0 )
        {
            mlpcreate1(nin, nhid1, nout, network);
        }
        else
        {
            if( nkind==1 )
            {
                mlpcreatec1(nin, nhid1, nout, network);
            }
            else
            {
                if( nkind==2 )
                {
                    mlpcreateb1(nin, nhid1, nout, a1, a2, network);
                }
                else
                {
                    if( nkind==3 )
                    {
                        mlpcreater1(nin, nhid1, nout, a1, a2, network);
                    }
                }
            }
        }
        return;
    }
    
    //
    // Two hidden layers
    //
    if( nkind==0 )
    {
        mlpcreate2(nin, nhid1, nhid2, nout, network);
    }
    else
    {
        if( nkind==1 )
        {
            mlpcreatec2(nin, nhid1, nhid2, nout, network);
        }
        else
        {
            if( nkind==2 )
            {
                mlpcreateb2(nin, nhid1, nhid2, nout, a1, a2, network);
            }
            else
            {
                if( nkind==3 )
                {
                    mlpcreater2(nin, nhid1, nhid2, nout, a1, a2, network);
                }
            }
        }
    }
}


/*************************************************************************
Unsets network (initialize it to smallest network possible
*************************************************************************/
static void unsetnetwork(multilayerperceptron& network)
{

    mlpcreate0(1, 1, network);
}


/*************************************************************************
Iformational functions test
*************************************************************************/
static void testinformational(int nkind,
     int nin,
     int nhid1,
     int nhid2,
     int nout,
     int passcount,
     bool& err)
{
    multilayerperceptron network;
    int n1;
    int n2;
    int wcount;

    createnetwork(network, nkind, 0.0, 0.0, nin, nhid1, nhid2, nout);
    mlpproperties(network, n1, n2, wcount);
    err = err||n1!=nin||n2!=nout||wcount<=0;
}


/*************************************************************************
Processing functions test
*************************************************************************/
static void testprocessing(int nkind,
     int nin,
     int nhid1,
     int nhid2,
     int nout,
     int passcount,
     bool& err)
{
    multilayerperceptron network;
    multilayerperceptron network2;
    int n1;
    int n2;
    int wcount;
    bool zeronet;
    double a1;
    double a2;
    int pass;
    int i;
    bool allsame;
    int rlen;
    ap::real_1d_array x1;
    ap::real_1d_array x2;
    ap::real_1d_array y1;
    ap::real_1d_array y2;
    ap::real_1d_array ra;
    ap::real_1d_array ra2;
    double v;

    ap::ap_error::make_assertion(passcount>=2, "PassCount<2!");
    
    //
    // Prepare network
    //
    a1 = 0;
    a2 = 0;
    if( nkind==2 )
    {
        a1 = 1000*ap::randomreal()-500;
        a2 = 2*ap::randomreal()-1;
    }
    if( nkind==3 )
    {
        a1 = 1000*ap::randomreal()-500;
        a2 = a1+(2*ap::randominteger(2)-1)*(0.1+0.9*ap::randomreal());
    }
    createnetwork(network, nkind, a1, a2, nin, nhid1, nhid2, nout);
    mlpproperties(network, n1, n2, wcount);
    
    //
    // Initialize arrays
    //
    x1.setbounds(0, nin-1);
    x2.setbounds(0, nin-1);
    y1.setbounds(0, nout-1);
    y2.setbounds(0, nout-1);
    
    //
    // Main cycle
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        
        //
        // Last run is made on zero network
        //
        mlprandomizefull(network);
        zeronet = false;
        if( pass==passcount )
        {
            ap::vmul(&network.weights(0), 1, ap::vlen(0,wcount-1), 0);
            zeronet = true;
        }
        
        //
        // Same inputs leads to same outputs
        //
        for(i = 0; i <= nin-1; i++)
        {
            x1(i) = 2*ap::randomreal()-1;
            x2(i) = x1(i);
        }
        for(i = 0; i <= nout-1; i++)
        {
            y1(i) = 2*ap::randomreal()-1;
            y2(i) = 2*ap::randomreal()-1;
        }
        mlpprocess(network, x1, y1);
        mlpprocess(network, x2, y2);
        allsame = true;
        for(i = 0; i <= nout-1; i++)
        {
            allsame = allsame&&ap::fp_eq(y1(i),y2(i));
        }
        err = err||!allsame;
        
        //
        // Same inputs on original network leads to same outputs
        // on copy created using MLPCopy
        //
        unsetnetwork(network2);
        mlpcopy(network, network2);
        for(i = 0; i <= nin-1; i++)
        {
            x1(i) = 2*ap::randomreal()-1;
            x2(i) = x1(i);
        }
        for(i = 0; i <= nout-1; i++)
        {
            y1(i) = 2*ap::randomreal()-1;
            y2(i) = 2*ap::randomreal()-1;
        }
        mlpprocess(network, x1, y1);
        mlpprocess(network2, x2, y2);
        allsame = true;
        for(i = 0; i <= nout-1; i++)
        {
            allsame = allsame&&ap::fp_eq(y1(i),y2(i));
        }
        err = err||!allsame;
        
        //
        // Same inputs on original network leads to same outputs
        // on copy created using MLPSerialize
        //
        unsetnetwork(network2);
        mlpserialize(network, ra, rlen);
        ra2.setbounds(0, rlen-1);
        for(i = 0; i <= rlen-1; i++)
        {
            ra2(i) = ra(i);
        }
        mlpunserialize(ra2, network2);
        for(i = 0; i <= nin-1; i++)
        {
            x1(i) = 2*ap::randomreal()-1;
            x2(i) = x1(i);
        }
        for(i = 0; i <= nout-1; i++)
        {
            y1(i) = 2*ap::randomreal()-1;
            y2(i) = 2*ap::randomreal()-1;
        }
        mlpprocess(network, x1, y1);
        mlpprocess(network2, x2, y2);
        allsame = true;
        for(i = 0; i <= nout-1; i++)
        {
            allsame = allsame&&ap::fp_eq(y1(i),y2(i));
        }
        err = err||!allsame;
        
        //
        // Different inputs leads to different outputs (non-zero network)
        //
        if( !zeronet )
        {
            for(i = 0; i <= nin-1; i++)
            {
                x1(i) = 2*ap::randomreal()-1;
                x2(i) = 2*ap::randomreal()-1;
            }
            for(i = 0; i <= nout-1; i++)
            {
                y1(i) = 2*ap::randomreal()-1;
                y2(i) = y1(i);
            }
            mlpprocess(network, x1, y1);
            mlpprocess(network, x2, y2);
            allsame = true;
            for(i = 0; i <= nout-1; i++)
            {
                allsame = allsame&&ap::fp_eq(y1(i),y2(i));
            }
            err = err||allsame;
        }
        
        //
        // Randomization changes outputs (when inputs are unchanged, non-zero network)
        //
        if( !zeronet )
        {
            for(i = 0; i <= nin-1; i++)
            {
                x1(i) = 2*ap::randomreal()-1;
                x2(i) = 2*ap::randomreal()-1;
            }
            for(i = 0; i <= nout-1; i++)
            {
                y1(i) = 2*ap::randomreal()-1;
                y2(i) = y1(i);
            }
            mlpcopy(network, network2);
            mlprandomize(network2);
            mlpprocess(network, x1, y1);
            mlpprocess(network2, x1, y2);
            allsame = true;
            for(i = 0; i <= nout-1; i++)
            {
                allsame = allsame&&ap::fp_eq(y1(i),y2(i));
            }
            err = err||allsame;
        }
        
        //
        // Full randomization changes outputs (when inputs are unchanged, non-zero network)
        //
        if( !zeronet )
        {
            for(i = 0; i <= nin-1; i++)
            {
                x1(i) = 2*ap::randomreal()-1;
                x2(i) = 2*ap::randomreal()-1;
            }
            for(i = 0; i <= nout-1; i++)
            {
                y1(i) = 2*ap::randomreal()-1;
                y2(i) = y1(i);
            }
            mlpcopy(network, network2);
            mlprandomizefull(network2);
            mlpprocess(network, x1, y1);
            mlpprocess(network2, x1, y2);
            allsame = true;
            for(i = 0; i <= nout-1; i++)
            {
                allsame = allsame&&ap::fp_eq(y1(i),y2(i));
            }
            err = err||allsame;
        }
        
        //
        // Normalization properties
        //
        if( nkind==1 )
        {
            
            //
            // Classifier network outputs are normalized
            //
            for(i = 0; i <= nin-1; i++)
            {
                x1(i) = 2*ap::randomreal()-1;
            }
            mlpprocess(network, x1, y1);
            v = 0;
            for(i = 0; i <= nout-1; i++)
            {
                v = v+y1(i);
                err = err||ap::fp_less(y1(i),0);
            }
            err = err||ap::fp_greater(fabs(v-1),1000*ap::machineepsilon);
        }
        if( nkind==2 )
        {
            
            //
            // B-type network outputs are bounded from above/below
            //
            for(i = 0; i <= nin-1; i++)
            {
                x1(i) = 2*ap::randomreal()-1;
            }
            mlpprocess(network, x1, y1);
            for(i = 0; i <= nout-1; i++)
            {
                if( ap::fp_greater_eq(a2,0) )
                {
                    err = err||ap::fp_less(y1(i),a1);
                }
                else
                {
                    err = err||ap::fp_greater(y1(i),a1);
                }
            }
        }
        if( nkind==3 )
        {
            
            //
            // R-type network outputs are within [A1,A2] (or [A2,A1])
            //
            for(i = 0; i <= nin-1; i++)
            {
                x1(i) = 2*ap::randomreal()-1;
            }
            mlpprocess(network, x1, y1);
            for(i = 0; i <= nout-1; i++)
            {
                err = err||ap::fp_less(y1(i),ap::minreal(a1, a2))||ap::fp_greater(y1(i),ap::maxreal(a1, a2));
            }
        }
    }
}


/*************************************************************************
Gradient functions test
*************************************************************************/
static void testgradient(int nkind,
     int nin,
     int nhid1,
     int nhid2,
     int nout,
     int passcount,
     bool& err)
{
    multilayerperceptron network;
    multilayerperceptron network2;
    int n1;
    int n2;
    int wcount;
    bool zeronet;
    double h;
    double etol;
    double a1;
    double a2;
    int pass;
    int i;
    int j;
    int k;
    bool allsame;
    int ilen;
    int rlen;
    int ssize;
    ap::real_2d_array xy;
    ap::real_1d_array grad1;
    ap::real_1d_array grad2;
    ap::real_1d_array x;
    ap::real_1d_array y;
    ap::real_1d_array x1;
    ap::real_1d_array x2;
    ap::real_1d_array y1;
    ap::real_1d_array y2;
    ap::integer_1d_array ia;
    ap::real_1d_array ra;
    double v;
    double e;
    double e1;
    double e2;
    double v1;
    double v2;
    double v3;
    double v4;
    double wprev;

    ap::ap_error::make_assertion(passcount>=2, "PassCount<2!");
    a1 = 0;
    a2 = 0;
    if( nkind==2 )
    {
        a1 = 1000*ap::randomreal()-500;
        a2 = 2*ap::randomreal()-1;
    }
    if( nkind==3 )
    {
        a1 = 1000*ap::randomreal()-500;
        a2 = a1+(2*ap::randominteger(2)-1)*(0.1+0.9*ap::randomreal());
    }
    createnetwork(network, nkind, a1, a2, nin, nhid1, nhid2, nout);
    mlpproperties(network, n1, n2, wcount);
    h = 0.0001;
    etol = 0.01;
    
    //
    // Initialize
    //
    x.setbounds(0, nin-1);
    x1.setbounds(0, nin-1);
    x2.setbounds(0, nin-1);
    y.setbounds(0, nout-1);
    y1.setbounds(0, nout-1);
    y2.setbounds(0, nout-1);
    grad1.setbounds(0, wcount-1);
    grad2.setbounds(0, wcount-1);
    
    //
    // Process
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        mlprandomizefull(network);
        
        //
        // Test error/gradient calculation (least squares)
        //
        xy.setbounds(0, 0, 0, nin+nout-1);
        for(i = 0; i <= nin-1; i++)
        {
            x(i) = 4*ap::randomreal()-2;
        }
        ap::vmove(&xy(0, 0), 1, &x(0), 1, ap::vlen(0,nin-1));
        if( mlpissoftmax(network) )
        {
            for(i = 0; i <= nout-1; i++)
            {
                y(i) = 0;
            }
            xy(0,nin) = ap::randominteger(nout);
            y(ap::round(xy(0,nin))) = 1;
        }
        else
        {
            for(i = 0; i <= nout-1; i++)
            {
                y(i) = 4*ap::randomreal()-2;
            }
            ap::vmove(&xy(0, nin), 1, &y(0), 1, ap::vlen(nin,nin+nout-1));
        }
        mlpgrad(network, x, y, e, grad2);
        mlpprocess(network, x, y2);
        ap::vsub(&y2(0), 1, &y(0), 1, ap::vlen(0,nout-1));
        v = ap::vdotproduct(&y2(0), 1, &y2(0), 1, ap::vlen(0,nout-1));
        v = v/2;
        err = err||ap::fp_greater(fabs((v-e)/v),etol);
        err = err||ap::fp_greater(fabs((mlperror(network, xy, 1)-v)/v),etol);
        for(i = 0; i <= wcount-1; i++)
        {
            wprev = network.weights(i);
            network.weights(i) = wprev-2*h;
            mlpprocess(network, x, y1);
            ap::vsub(&y1(0), 1, &y(0), 1, ap::vlen(0,nout-1));
            v1 = ap::vdotproduct(&y1(0), 1, &y1(0), 1, ap::vlen(0,nout-1));
            v1 = v1/2;
            network.weights(i) = wprev-h;
            mlpprocess(network, x, y1);
            ap::vsub(&y1(0), 1, &y(0), 1, ap::vlen(0,nout-1));
            v2 = ap::vdotproduct(&y1(0), 1, &y1(0), 1, ap::vlen(0,nout-1));
            v2 = v2/2;
            network.weights(i) = wprev+h;
            mlpprocess(network, x, y1);
            ap::vsub(&y1(0), 1, &y(0), 1, ap::vlen(0,nout-1));
            v3 = ap::vdotproduct(&y1(0), 1, &y1(0), 1, ap::vlen(0,nout-1));
            v3 = v3/2;
            network.weights(i) = wprev+2*h;
            mlpprocess(network, x, y1);
            ap::vsub(&y1(0), 1, &y(0), 1, ap::vlen(0,nout-1));
            v4 = ap::vdotproduct(&y1(0), 1, &y1(0), 1, ap::vlen(0,nout-1));
            v4 = v4/2;
            network.weights(i) = wprev;
            grad1(i) = (v1-8*v2+8*v3-v4)/(12*h);
            if( ap::fp_greater(fabs(grad1(i)),1.0E-3) )
            {
                err = err||ap::fp_greater(fabs((grad2(i)-grad1(i))/grad1(i)),etol);
            }
            else
            {
                err = err||ap::fp_greater(fabs(grad2(i)-grad1(i)),etol);
            }
        }
        
        //
        // Test error/gradient calculation (natural).
        // Testing on non-random structure networks
        // (because NKind is representative only in that case).
        //
        xy.setbounds(0, 0, 0, nin+nout-1);
        for(i = 0; i <= nin-1; i++)
        {
            x(i) = 4*ap::randomreal()-2;
        }
        ap::vmove(&xy(0, 0), 1, &x(0), 1, ap::vlen(0,nin-1));
        if( mlpissoftmax(network) )
        {
            for(i = 0; i <= nout-1; i++)
            {
                y(i) = 0;
            }
            xy(0,nin) = ap::randominteger(nout);
            y(ap::round(xy(0,nin))) = 1;
        }
        else
        {
            for(i = 0; i <= nout-1; i++)
            {
                y(i) = 4*ap::randomreal()-2;
            }
            ap::vmove(&xy(0, nin), 1, &y(0), 1, ap::vlen(nin,nin+nout-1));
        }
        mlpgradn(network, x, y, e, grad2);
        mlpprocess(network, x, y2);
        v = 0;
        if( nkind!=1 )
        {
            for(i = 0; i <= nout-1; i++)
            {
                v = v+0.5*ap::sqr(y2(i)-y(i));
            }
        }
        else
        {
            for(i = 0; i <= nout-1; i++)
            {
                if( ap::fp_neq(y(i),0) )
                {
                    if( ap::fp_eq(y2(i),0) )
                    {
                        v = v+y(i)*log(ap::maxrealnumber);
                    }
                    else
                    {
                        v = v+y(i)*log(y(i)/y2(i));
                    }
                }
            }
        }
        err = err||ap::fp_greater(fabs((v-e)/v),etol);
        err = err||ap::fp_greater(fabs((mlperrorn(network, xy, 1)-v)/v),etol);
        for(i = 0; i <= wcount-1; i++)
        {
            wprev = network.weights(i);
            network.weights(i) = wprev+h;
            mlpprocess(network, x, y2);
            network.weights(i) = wprev-h;
            mlpprocess(network, x, y1);
            network.weights(i) = wprev;
            v = 0;
            if( nkind!=1 )
            {
                for(j = 0; j <= nout-1; j++)
                {
                    v = v+0.5*(ap::sqr(y2(j)-y(j))-ap::sqr(y1(j)-y(j)))/(2*h);
                }
            }
            else
            {
                for(j = 0; j <= nout-1; j++)
                {
                    if( ap::fp_neq(y(j),0) )
                    {
                        if( ap::fp_eq(y2(j),0) )
                        {
                            v = v+y(j)*log(ap::maxrealnumber);
                        }
                        else
                        {
                            v = v+y(j)*log(y(j)/y2(j));
                        }
                        if( ap::fp_eq(y1(j),0) )
                        {
                            v = v-y(j)*log(ap::maxrealnumber);
                        }
                        else
                        {
                            v = v-y(j)*log(y(j)/y1(j));
                        }
                    }
                }
                v = v/(2*h);
            }
            grad1(i) = v;
            if( ap::fp_greater(fabs(grad1(i)),1.0E-3) )
            {
                err = err||ap::fp_greater(fabs((grad2(i)-grad1(i))/grad1(i)),etol);
            }
            else
            {
                err = err||ap::fp_greater(fabs(grad2(i)-grad1(i)),etol);
            }
        }
        
        //
        // Test gradient calculation: batch (least squares)
        //
        ssize = 1+ap::randominteger(10);
        xy.setbounds(0, ssize-1, 0, nin+nout-1);
        for(i = 0; i <= wcount-1; i++)
        {
            grad1(i) = 0;
        }
        e1 = 0;
        for(i = 0; i <= ssize-1; i++)
        {
            for(j = 0; j <= nin-1; j++)
            {
                x1(j) = 4*ap::randomreal()-2;
            }
            ap::vmove(&xy(i, 0), 1, &x1(0), 1, ap::vlen(0,nin-1));
            if( mlpissoftmax(network) )
            {
                for(j = 0; j <= nout-1; j++)
                {
                    y1(j) = 0;
                }
                xy(i,nin) = ap::randominteger(nout);
                y1(ap::round(xy(i,nin))) = 1;
            }
            else
            {
                for(j = 0; j <= nout-1; j++)
                {
                    y1(j) = 4*ap::randomreal()-2;
                }
                ap::vmove(&xy(i, nin), 1, &y1(0), 1, ap::vlen(nin,nin+nout-1));
            }
            mlpgrad(network, x1, y1, v, grad2);
            e1 = e1+v;
            ap::vadd(&grad1(0), 1, &grad2(0), 1, ap::vlen(0,wcount-1));
        }
        mlpgradbatch(network, xy, ssize, e2, grad2);
        err = err||ap::fp_greater(fabs(e1-e2)/e1,0.01);
        for(i = 0; i <= wcount-1; i++)
        {
            if( ap::fp_neq(grad1(i),0) )
            {
                err = err||ap::fp_greater(fabs((grad2(i)-grad1(i))/grad1(i)),etol);
            }
            else
            {
                err = err||ap::fp_neq(grad2(i),grad1(i));
            }
        }
        
        //
        // Test gradient calculation: batch (natural error func)
        //
        ssize = 1+ap::randominteger(10);
        xy.setbounds(0, ssize-1, 0, nin+nout-1);
        for(i = 0; i <= wcount-1; i++)
        {
            grad1(i) = 0;
        }
        e1 = 0;
        for(i = 0; i <= ssize-1; i++)
        {
            for(j = 0; j <= nin-1; j++)
            {
                x1(j) = 4*ap::randomreal()-2;
            }
            ap::vmove(&xy(i, 0), 1, &x1(0), 1, ap::vlen(0,nin-1));
            if( mlpissoftmax(network) )
            {
                for(j = 0; j <= nout-1; j++)
                {
                    y1(j) = 0;
                }
                xy(i,nin) = ap::randominteger(nout);
                y1(ap::round(xy(i,nin))) = 1;
            }
            else
            {
                for(j = 0; j <= nout-1; j++)
                {
                    y1(j) = 4*ap::randomreal()-2;
                }
                ap::vmove(&xy(i, nin), 1, &y1(0), 1, ap::vlen(nin,nin+nout-1));
            }
            mlpgradn(network, x1, y1, v, grad2);
            e1 = e1+v;
            ap::vadd(&grad1(0), 1, &grad2(0), 1, ap::vlen(0,wcount-1));
        }
        mlpgradnbatch(network, xy, ssize, e2, grad2);
        err = err||ap::fp_greater(fabs(e1-e2)/e1,etol);
        for(i = 0; i <= wcount-1; i++)
        {
            if( ap::fp_neq(grad1(i),0) )
            {
                err = err||ap::fp_greater(fabs((grad2(i)-grad1(i))/grad1(i)),etol);
            }
            else
            {
                err = err||ap::fp_neq(grad2(i),grad1(i));
            }
        }
    }
}


/*************************************************************************
Hessian functions test
*************************************************************************/
static void testhessian(int nkind,
     int nin,
     int nhid1,
     int nhid2,
     int nout,
     int passcount,
     bool& err)
{
    multilayerperceptron network;
    multilayerperceptron network2;
    int hkind;
    int n1;
    int n2;
    int wcount;
    bool zeronet;
    double h;
    double etol;
    int pass;
    int i;
    int j;
    bool allsame;
    int ilen;
    int rlen;
    int ssize;
    double a1;
    double a2;
    ap::real_2d_array xy;
    ap::real_2d_array h1;
    ap::real_2d_array h2;
    ap::real_1d_array grad1;
    ap::real_1d_array grad2;
    ap::real_1d_array grad3;
    ap::real_1d_array x;
    ap::real_1d_array y;
    ap::real_1d_array x1;
    ap::real_1d_array x2;
    ap::real_1d_array y1;
    ap::real_1d_array y2;
    ap::integer_1d_array ia;
    ap::real_1d_array ra;
    double v;
    double e;
    double e1;
    double e2;
    double v1;
    double v2;
    double v3;
    double v4;
    double wprev;

    ap::ap_error::make_assertion(passcount>=2, "PassCount<2!");
    a1 = 0;
    a2 = 0;
    if( nkind==2 )
    {
        a1 = 1000*ap::randomreal()-500;
        a2 = 2*ap::randomreal()-1;
    }
    if( nkind==3 )
    {
        a1 = 1000*ap::randomreal()-500;
        a2 = a1+(2*ap::randominteger(2)-1)*(0.1+0.9*ap::randomreal());
    }
    createnetwork(network, nkind, a1, a2, nin, nhid1, nhid2, nout);
    mlpproperties(network, n1, n2, wcount);
    h = 0.0001;
    etol = 0.05;
    
    //
    // Initialize
    //
    x.setbounds(0, nin-1);
    x1.setbounds(0, nin-1);
    x2.setbounds(0, nin-1);
    y.setbounds(0, nout-1);
    y1.setbounds(0, nout-1);
    y2.setbounds(0, nout-1);
    grad1.setbounds(0, wcount-1);
    grad2.setbounds(0, wcount-1);
    grad3.setbounds(0, wcount-1);
    h1.setbounds(0, wcount-1, 0, wcount-1);
    h2.setbounds(0, wcount-1, 0, wcount-1);
    
    //
    // Process
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        mlprandomizefull(network);
        
        //
        // Test hessian calculation .
        // E1 contains total error (calculated using MLPGrad/MLPGradN)
        // Grad1 contains total gradient (calculated using MLPGrad/MLPGradN)
        // H1 contains Hessian calculated using differences of gradients
        //
        // E2, Grad2 and H2 contains corresponing values calculated using MLPHessianBatch/MLPHessianNBatch
        //
        for(hkind = 0; hkind <= 1; hkind++)
        {
            ssize = 1+ap::randominteger(10);
            xy.setbounds(0, ssize-1, 0, nin+nout-1);
            for(i = 0; i <= wcount-1; i++)
            {
                grad1(i) = 0;
            }
            for(i = 0; i <= wcount-1; i++)
            {
                for(j = 0; j <= wcount-1; j++)
                {
                    h1(i,j) = 0;
                }
            }
            e1 = 0;
            for(i = 0; i <= ssize-1; i++)
            {
                
                //
                // X, Y
                //
                for(j = 0; j <= nin-1; j++)
                {
                    x1(j) = 4*ap::randomreal()-2;
                }
                ap::vmove(&xy(i, 0), 1, &x1(0), 1, ap::vlen(0,nin-1));
                if( mlpissoftmax(network) )
                {
                    for(j = 0; j <= nout-1; j++)
                    {
                        y1(j) = 0;
                    }
                    xy(i,nin) = ap::randominteger(nout);
                    y1(ap::round(xy(i,nin))) = 1;
                }
                else
                {
                    for(j = 0; j <= nout-1; j++)
                    {
                        y1(j) = 4*ap::randomreal()-2;
                    }
                    ap::vmove(&xy(i, nin), 1, &y1(0), 1, ap::vlen(nin,nin+nout-1));
                }
                
                //
                // E1, Grad1
                //
                if( hkind==0 )
                {
                    mlpgrad(network, x1, y1, v, grad2);
                }
                else
                {
                    mlpgradn(network, x1, y1, v, grad2);
                }
                e1 = e1+v;
                ap::vadd(&grad1(0), 1, &grad2(0), 1, ap::vlen(0,wcount-1));
                
                //
                // H1
                //
                for(j = 0; j <= wcount-1; j++)
                {
                    wprev = network.weights(j);
                    network.weights(j) = wprev-2*h;
                    if( hkind==0 )
                    {
                        mlpgrad(network, x1, y1, v, grad2);
                    }
                    else
                    {
                        mlpgradn(network, x1, y1, v, grad2);
                    }
                    network.weights(j) = wprev-h;
                    if( hkind==0 )
                    {
                        mlpgrad(network, x1, y1, v, grad3);
                    }
                    else
                    {
                        mlpgradn(network, x1, y1, v, grad3);
                    }
                    ap::vsub(&grad2(0), 1, &grad3(0), 1, ap::vlen(0,wcount-1), 8);
                    network.weights(j) = wprev+h;
                    if( hkind==0 )
                    {
                        mlpgrad(network, x1, y1, v, grad3);
                    }
                    else
                    {
                        mlpgradn(network, x1, y1, v, grad3);
                    }
                    ap::vadd(&grad2(0), 1, &grad3(0), 1, ap::vlen(0,wcount-1), 8);
                    network.weights(j) = wprev+2*h;
                    if( hkind==0 )
                    {
                        mlpgrad(network, x1, y1, v, grad3);
                    }
                    else
                    {
                        mlpgradn(network, x1, y1, v, grad3);
                    }
                    ap::vsub(&grad2(0), 1, &grad3(0), 1, ap::vlen(0,wcount-1));
                    v = 1/(12*h);
                    ap::vadd(&h1(j, 0), 1, &grad2(0), 1, ap::vlen(0,wcount-1), v);
                    network.weights(j) = wprev;
                }
            }
            if( hkind==0 )
            {
                mlphessianbatch(network, xy, ssize, e2, grad2, h2);
            }
            else
            {
                mlphessiannbatch(network, xy, ssize, e2, grad2, h2);
            }
            err = err||ap::fp_greater(fabs(e1-e2)/e1,etol);
            for(i = 0; i <= wcount-1; i++)
            {
                if( ap::fp_greater(fabs(grad1(i)),1.0E-2) )
                {
                    err = err||ap::fp_greater(fabs((grad2(i)-grad1(i))/grad1(i)),etol);
                }
                else
                {
                    err = err||ap::fp_greater(fabs(grad2(i)-grad1(i)),etol);
                }
            }
            for(i = 0; i <= wcount-1; i++)
            {
                for(j = 0; j <= wcount-1; j++)
                {
                    if( ap::fp_greater(fabs(h1(i,j)),5.0E-2) )
                    {
                        err = err||ap::fp_greater(fabs((h1(i,j)-h2(i,j))/h1(i,j)),etol);
                    }
                    else
                    {
                        err = err||ap::fp_greater(fabs(h2(i,j)-h1(i,j)),etol);
                    }
                }
            }
        }
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testmlpunit_test_silent()
{
    bool result;

    result = testmlp(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testmlpunit_test()
{
    bool result;

    result = testmlp(false);
    return result;
}





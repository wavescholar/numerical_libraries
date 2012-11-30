
 
#include <stdio.h>
#include "testmlpeunit.h"

static void createensemble(mlpensemble& ensemble,
     int nkind,
     double a1,
     double a2,
     int nin,
     int nhid1,
     int nhid2,
     int nout,
     int ec);
static void unsetensemble(mlpensemble& ensemble);
static void testinformational(int nkind,
     int nin,
     int nhid1,
     int nhid2,
     int nout,
     int ec,
     int passcount,
     bool& err);
static void testprocessing(int nkind,
     int nin,
     int nhid1,
     int nhid2,
     int nout,
     int ec,
     int passcount,
     bool& err);

bool testmlpe(bool silent)
{
    bool result;
    bool waserrors;
    int passcount;
    int maxn;
    int maxhid;
    int nf;
    int nhid;
    int nl;
    int nhid1;
    int nhid2;
    int ec;
    int nkind;
    int algtype;
    int tasktype;
    int pass;
    mlpensemble ensemble;
    mlpreport rep;
    mlpcvreport oobrep;
    ap::real_2d_array xy;
    int i;
    int j;
    int nin;
    int nout;
    int npoints;
    double e;
    int info;
    int nless;
    int nall;
    int nclasses;
    bool allsame;
    bool inferrors;
    bool procerrors;
    bool trnerrors;

    waserrors = false;
    inferrors = false;
    procerrors = false;
    trnerrors = false;
    passcount = 10;
    maxn = 4;
    maxhid = 4;
    
    //
    // General MLP ensembles tests
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
                        for(ec = 1; ec <= 3; ec++)
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
                            testinformational(nkind, nf, nhid1, nhid2, nl, ec, passcount, inferrors);
                            testprocessing(nkind, nf, nhid1, nhid2, nl, ec, passcount, procerrors);
                        }
                    }
                }
            }
        }
    }
    
    //
    // network training must reduce error
    // test on random regression task
    //
    nin = 3;
    nout = 2;
    nhid = 5;
    npoints = 100;
    nless = 0;
    nall = 0;
    for(pass = 1; pass <= 10; pass++)
    {
        for(algtype = 0; algtype <= 1; algtype++)
        {
            for(tasktype = 0; tasktype <= 1; tasktype++)
            {
                if( tasktype==0 )
                {
                    xy.setbounds(0, npoints-1, 0, nin+nout-1);
                    for(i = 0; i <= npoints-1; i++)
                    {
                        for(j = 0; j <= nin+nout-1; j++)
                        {
                            xy(i,j) = 2*ap::randomreal()-1;
                        }
                    }
                    mlpecreate1(nin, nhid, nout, 1+ap::randominteger(3), ensemble);
                }
                else
                {
                    xy.setbounds(0, npoints-1, 0, nin);
                    nclasses = 2+ap::randominteger(2);
                    for(i = 0; i <= npoints-1; i++)
                    {
                        for(j = 0; j <= nin-1; j++)
                        {
                            xy(i,j) = 2*ap::randomreal()-1;
                        }
                        xy(i,nin) = ap::randominteger(nclasses);
                    }
                    mlpecreatec1(nin, nhid, nclasses, 1+ap::randominteger(3), ensemble);
                }
                e = mlpermserror(ensemble, xy, npoints);
                if( algtype==0 )
                {
                    mlpebagginglm(ensemble, xy, npoints, 0.001, 1, info, rep, oobrep);
                }
                else
                {
                    mlpebagginglbfgs(ensemble, xy, npoints, 0.001, 1, 0.01, 0, info, rep, oobrep);
                }
                if( info<0 )
                {
                    trnerrors = true;
                }
                else
                {
                    if( ap::fp_less(mlpermserror(ensemble, xy, npoints),e) )
                    {
                        nless = nless+1;
                    }
                }
                nall = nall+1;
            }
        }
    }
    trnerrors = trnerrors||ap::fp_greater(nall-nless,0.3*nall);
    
    //
    // Final report
    //
    waserrors = inferrors||procerrors||trnerrors;
    if( !silent )
    {
        printf("MLP ENSEMBLE TEST\n");
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
static void createensemble(mlpensemble& ensemble,
     int nkind,
     double a1,
     double a2,
     int nin,
     int nhid1,
     int nhid2,
     int nout,
     int ec)
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
            mlpecreate0(nin, nout, ec, ensemble);
        }
        else
        {
            if( nkind==1 )
            {
                mlpecreatec0(nin, nout, ec, ensemble);
            }
            else
            {
                if( nkind==2 )
                {
                    mlpecreateb0(nin, nout, a1, a2, ec, ensemble);
                }
                else
                {
                    if( nkind==3 )
                    {
                        mlpecreater0(nin, nout, a1, a2, ec, ensemble);
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
            mlpecreate1(nin, nhid1, nout, ec, ensemble);
        }
        else
        {
            if( nkind==1 )
            {
                mlpecreatec1(nin, nhid1, nout, ec, ensemble);
            }
            else
            {
                if( nkind==2 )
                {
                    mlpecreateb1(nin, nhid1, nout, a1, a2, ec, ensemble);
                }
                else
                {
                    if( nkind==3 )
                    {
                        mlpecreater1(nin, nhid1, nout, a1, a2, ec, ensemble);
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
        mlpecreate2(nin, nhid1, nhid2, nout, ec, ensemble);
    }
    else
    {
        if( nkind==1 )
        {
            mlpecreatec2(nin, nhid1, nhid2, nout, ec, ensemble);
        }
        else
        {
            if( nkind==2 )
            {
                mlpecreateb2(nin, nhid1, nhid2, nout, a1, a2, ec, ensemble);
            }
            else
            {
                if( nkind==3 )
                {
                    mlpecreater2(nin, nhid1, nhid2, nout, a1, a2, ec, ensemble);
                }
            }
        }
    }
}


/*************************************************************************
Unsets network (initialize it to smallest network possible
*************************************************************************/
static void unsetensemble(mlpensemble& ensemble)
{

    mlpecreate0(1, 1, 1, ensemble);
}


/*************************************************************************
Iformational functions test
*************************************************************************/
static void testinformational(int nkind,
     int nin,
     int nhid1,
     int nhid2,
     int nout,
     int ec,
     int passcount,
     bool& err)
{
    mlpensemble ensemble;
    int n1;
    int n2;

    createensemble(ensemble, nkind, -1.0, 1.0, nin, nhid1, nhid2, nout, ec);
    mlpeproperties(ensemble, n1, n2);
    err = err||n1!=nin||n2!=nout;
}


/*************************************************************************
Processing functions test
*************************************************************************/
static void testprocessing(int nkind,
     int nin,
     int nhid1,
     int nhid2,
     int nout,
     int ec,
     int passcount,
     bool& err)
{
    mlpensemble ensemble;
    mlpensemble ensemble2;
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
        createensemble(ensemble, nkind, a1, a2, nin, nhid1, nhid2, nout, ec);
        
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
        mlpeprocess(ensemble, x1, y1);
        mlpeprocess(ensemble, x2, y2);
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
        unsetensemble(ensemble2);
        mlpecopy(ensemble, ensemble2);
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
        mlpeprocess(ensemble, x1, y1);
        mlpeprocess(ensemble2, x2, y2);
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
        unsetensemble(ensemble2);
        mlpeserialize(ensemble, ra, rlen);
        ra2.setbounds(0, rlen-1);
        for(i = 0; i <= rlen-1; i++)
        {
            ra2(i) = ra(i);
        }
        mlpeunserialize(ra2, ensemble2);
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
        mlpeprocess(ensemble, x1, y1);
        mlpeprocess(ensemble2, x2, y2);
        allsame = true;
        for(i = 0; i <= nout-1; i++)
        {
            allsame = allsame&&ap::fp_eq(y1(i),y2(i));
        }
        err = err||!allsame;
        
        //
        // Different inputs leads to different outputs (non-zero network)
        //
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
        mlpeprocess(ensemble, x1, y1);
        mlpeprocess(ensemble, x2, y2);
        allsame = true;
        for(i = 0; i <= nout-1; i++)
        {
            allsame = allsame&&ap::fp_eq(y1(i),y2(i));
        }
        err = err||allsame;
        
        //
        // Randomization changes outputs (when inputs are unchanged, non-zero network)
        //
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
        mlpecopy(ensemble, ensemble2);
        mlperandomize(ensemble2);
        mlpeprocess(ensemble, x1, y1);
        mlpeprocess(ensemble2, x1, y2);
        allsame = true;
        for(i = 0; i <= nout-1; i++)
        {
            allsame = allsame&&ap::fp_eq(y1(i),y2(i));
        }
        err = err||allsame;
        
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
            mlpeprocess(ensemble, x1, y1);
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
            mlpeprocess(ensemble, x1, y1);
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
            mlpeprocess(ensemble, x1, y1);
            for(i = 0; i <= nout-1; i++)
            {
                err = err||ap::fp_less(y1(i),ap::minreal(a1, a2))||ap::fp_greater(y1(i),ap::maxreal(a1, a2));
            }
        }
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testmlpeunit_test_silent()
{
    bool result;

    result = testmlpe(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testmlpeunit_test()
{
    bool result;

    result = testmlpe(false);
    return result;
}





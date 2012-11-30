
 
#include <stdio.h>
#include "testnearestneighborunit.h"

static void unset2d(ap::complex_2d_array& a);
static void unset1d(ap::real_1d_array& a);
static bool kdtresultsdifferent(const ap::real_2d_array& refxy,
     int ntotal,
     const ap::real_2d_array& qx,
     const ap::real_2d_array& qxy,
     const ap::integer_1d_array& qt,
     int n,
     int nx,
     int ny);
static double vnorm(const ap::real_1d_array& x, int n, int normtype);
static void testkdtuniform(const ap::real_2d_array& xy,
     const int& n,
     const int& nx,
     const int& ny,
     const int& normtype,
     bool& kdterrors);

/*************************************************************************
Testing Nearest Neighbor Search
*************************************************************************/
bool testnearestneighbor(bool silent)
{
    bool result;
    ap::real_2d_array xy;
    int i;
    int j;
    double v;
    int normtype;
    int nx;
    int ny;
    int n;
    int smalln;
    int largen;
    int passcount;
    int pass;
    bool waserrors;
    bool kdterrors;

    kdterrors = false;
    passcount = 2;
    smalln = 256;
    largen = 2048;
    ny = 3;
    
    //
    //
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        for(normtype = 0; normtype <= 2; normtype++)
        {
            for(nx = 1; nx <= 3; nx++)
            {
                
                //
                // Test in hypercube
                //
                xy.setlength(largen, nx+ny);
                for(i = 0; i <= largen-1; i++)
                {
                    for(j = 0; j <= nx+ny-1; j++)
                    {
                        xy(i,j) = 10*ap::randomreal()-5;
                    }
                }
                for(n = 1; n <= 10; n++)
                {
                    testkdtuniform(xy, n, nx, ap::randominteger(ny+1), normtype, kdterrors);
                }
                testkdtuniform(xy, largen, nx, ap::randominteger(ny+1), normtype, kdterrors);
                
                //
                // Test clustered (2*N points, pairs of equal points)
                //
                xy.setlength(2*smalln, nx+ny);
                for(i = 0; i <= smalln-1; i++)
                {
                    for(j = 0; j <= nx+ny-1; j++)
                    {
                        xy(2*i+0,j) = 10*ap::randomreal()-5;
                        xy(2*i+1,j) = xy(2*i+0,j);
                    }
                }
                testkdtuniform(xy, 2*smalln, nx, ap::randominteger(ny+1), normtype, kdterrors);
                
                //
                // Test degenerate case: all points are same except for one
                //
                xy.setlength(smalln, nx+ny);
                v = ap::randomreal();
                for(i = 0; i <= smalln-2; i++)
                {
                    for(j = 0; j <= nx+ny-1; j++)
                    {
                        xy(i,j) = v;
                    }
                }
                for(j = 0; j <= nx+ny-1; j++)
                {
                    xy(smalln-1,j) = 10*ap::randomreal()-5;
                }
                testkdtuniform(xy, smalln, nx, ap::randominteger(ny+1), normtype, kdterrors);
            }
        }
    }
    
    //
    // report
    //
    waserrors = kdterrors;
    if( !silent )
    {
        printf("TESTING NEAREST NEIGHBOR SEARCH\n");
        printf("* KD TREES:                              ");
        if( !kdterrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        if( waserrors )
        {
            printf("TEST FAILED\n");
        }
        else
        {
            printf("TEST PASSED\n");
        }
        printf("\n\n");
    }
    result = !waserrors;
    return result;
}


/*************************************************************************
Unsets 2D array.
*************************************************************************/
static void unset2d(ap::complex_2d_array& a)
{

    a.setbounds(0, 0, 0, 0);
    a(0,0) = 2*ap::randomreal()-1;
}


/*************************************************************************
Unsets 1D array.
*************************************************************************/
static void unset1d(ap::real_1d_array& a)
{

    a.setbounds(0, 0);
    a(0) = 2*ap::randomreal()-1;
}


/*************************************************************************
Compare results from different queries:
* X     just X-values
* XY    X-values and Y-values
* XT    X-values and tag values
*************************************************************************/
static bool kdtresultsdifferent(const ap::real_2d_array& refxy,
     int ntotal,
     const ap::real_2d_array& qx,
     const ap::real_2d_array& qxy,
     const ap::integer_1d_array& qt,
     int n,
     int nx,
     int ny)
{
    bool result;
    int i;
    int j;

    result = false;
    for(i = 0; i <= n-1; i++)
    {
        if( qt(i)<0||qt(i)>=ntotal )
        {
            result = true;
            return result;
        }
        for(j = 0; j <= nx-1; j++)
        {
            result = result||ap::fp_neq(qx(i,j),refxy(qt(i),j));
            result = result||ap::fp_neq(qxy(i,j),refxy(qt(i),j));
        }
        for(j = 0; j <= ny-1; j++)
        {
            result = result||ap::fp_neq(qxy(i,nx+j),refxy(qt(i),nx+j));
        }
    }
    return result;
}


/*************************************************************************
Returns norm
*************************************************************************/
static double vnorm(const ap::real_1d_array& x, int n, int normtype)
{
    double result;
    int i;

    result = ap::randomreal();
    if( normtype==0 )
    {
        result = 0;
        for(i = 0; i <= n-1; i++)
        {
            result = ap::maxreal(result, fabs(x(i)));
        }
        return result;
    }
    if( normtype==1 )
    {
        result = 0;
        for(i = 0; i <= n-1; i++)
        {
            result = result+fabs(x(i));
        }
        return result;
    }
    if( normtype==2 )
    {
        result = 0;
        for(i = 0; i <= n-1; i++)
        {
            result = result+ap::sqr(x(i));
        }
        result = sqrt(result);
        return result;
    }
    return result;
}


/*************************************************************************
Testing Nearest Neighbor Search on uniformly distributed hypercube

NormType: 0, 1, 2
D: space dimension
N: points count
*************************************************************************/
static void testkdtuniform(const ap::real_2d_array& xy,
     const int& n,
     const int& nx,
     const int& ny,
     const int& normtype,
     bool& kdterrors)
{
    double errtol;
    ap::integer_1d_array tags;
    ap::real_1d_array ptx;
    ap::real_1d_array tmpx;
    ap::boolean_1d_array tmpb;
    kdtree treex;
    kdtree treexy;
    kdtree treext;
    ap::real_2d_array qx;
    ap::real_2d_array qxy;
    ap::integer_1d_array qtags;
    ap::real_1d_array qr;
    int kx;
    int kxy;
    int kt;
    int kr;
    double eps;
    int i;
    int j;
    int k;
    int task;
    bool isequal;
    double r;
    int q;
    int qcount;

    qcount = 10;
    
    //
    // Tol - roundoff error tolerance (for '>=' comparisons)
    //
    errtol = 100000*ap::machineepsilon;
    
    //
    // fill tags
    //
    tags.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        tags(i) = i;
    }
    
    //
    // build trees
    //
    kdtreebuild(xy, n, nx, 0, normtype, treex);
    kdtreebuild(xy, n, nx, ny, normtype, treexy);
    kdtreebuildtagged(xy, tags, n, nx, 0, normtype, treext);
    
    //
    // allocate arrays
    //
    tmpx.setlength(nx);
    tmpb.setlength(n);
    qx.setlength(n, nx);
    qxy.setlength(n, nx+ny);
    qtags.setlength(n);
    qr.setlength(n);
    ptx.setlength(nx);
    
    //
    // test general K-NN queries (with self-matches):
    // * compare results from different trees (must be equal) and
    //   check that correct (value,tag) pairs are returned
    // * test results from XT tree - let R be radius of query result.
    //   then all points not in result must be not closer than R.
    //
    for(q = 1; q <= qcount; q++)
    {
        
        //
        // Select K: 1..N
        //
        if( ap::fp_greater(ap::randomreal(),0.5) )
        {
            k = 1+ap::randominteger(n);
        }
        else
        {
            k = 1;
        }
        
        //
        // Select point (either one of the points, or random)
        //
        if( ap::fp_greater(ap::randomreal(),0.5) )
        {
            i = ap::randominteger(n);
            ap::vmove(&ptx(0), 1, &xy(i, 0), 1, ap::vlen(0,nx-1));
        }
        else
        {
            for(i = 0; i <= nx-1; i++)
            {
                ptx(i) = 2*ap::randomreal()-1;
            }
        }
        
        //
        // Test:
        // * consistency of results from different queries
        // * points in query are IN the R-sphere (or at the boundary),
        //   and points not in query are outside of the R-sphere (or at the boundary)
        // * distances are correct and are ordered
        //
        kx = kdtreequeryknn(treex, ptx, k, true);
        kxy = kdtreequeryknn(treexy, ptx, k, true);
        kt = kdtreequeryknn(treext, ptx, k, true);
        if( kx!=k||kxy!=k||kt!=k )
        {
            kdterrors = true;
            return;
        }
        kx = 0;
        kxy = 0;
        kt = 0;
        kdtreequeryresultsx(treex, qx, kx);
        kdtreequeryresultsxy(treexy, qxy, kxy);
        kdtreequeryresultstags(treext, qtags, kt);
        kdtreequeryresultsdistances(treext, qr, kr);
        if( kx!=k||kxy!=k||kt!=k||kr!=k )
        {
            kdterrors = true;
            return;
        }
        kdterrors = kdterrors||kdtresultsdifferent(xy, n, qx, qxy, qtags, k, nx, ny);
        for(i = 0; i <= n-1; i++)
        {
            tmpb(i) = true;
        }
        r = 0;
        for(i = 0; i <= k-1; i++)
        {
            tmpb(qtags(i)) = false;
            ap::vmove(&tmpx(0), 1, &ptx(0), 1, ap::vlen(0,nx-1));
            ap::vsub(&tmpx(0), 1, &qx(i, 0), 1, ap::vlen(0,nx-1));
            r = ap::maxreal(r, vnorm(tmpx, nx, normtype));
        }
        for(i = 0; i <= n-1; i++)
        {
            if( tmpb(i) )
            {
                ap::vmove(&tmpx(0), 1, &ptx(0), 1, ap::vlen(0,nx-1));
                ap::vsub(&tmpx(0), 1, &xy(i, 0), 1, ap::vlen(0,nx-1));
                kdterrors = kdterrors||ap::fp_less(vnorm(tmpx, nx, normtype),r*(1-errtol));
            }
        }
        for(i = 0; i <= k-2; i++)
        {
            kdterrors = kdterrors||ap::fp_greater(qr(i),qr(i+1));
        }
        for(i = 0; i <= k-1; i++)
        {
            ap::vmove(&tmpx(0), 1, &ptx(0), 1, ap::vlen(0,nx-1));
            ap::vsub(&tmpx(0), 1, &xy(qtags(i), 0), 1, ap::vlen(0,nx-1));
            kdterrors = kdterrors||ap::fp_greater(fabs(vnorm(tmpx, nx, normtype)-qr(i)),errtol);
        }
    }
    
    //
    // test general approximate K-NN queries (with self-matches):
    // * compare results from different trees (must be equal) and
    //   check that correct (value,tag) pairs are returned
    // * test results from XT tree - let R be radius of query result.
    //   then all points not in result must be not closer than R/(1+Eps).
    //
    for(q = 1; q <= qcount; q++)
    {
        
        //
        // Select K: 1..N
        //
        if( ap::fp_greater(ap::randomreal(),0.5) )
        {
            k = 1+ap::randominteger(n);
        }
        else
        {
            k = 1;
        }
        
        //
        // Select Eps
        //
        eps = 0.5+ap::randomreal();
        
        //
        // Select point (either one of the points, or random)
        //
        if( ap::fp_greater(ap::randomreal(),0.5) )
        {
            i = ap::randominteger(n);
            ap::vmove(&ptx(0), 1, &xy(i, 0), 1, ap::vlen(0,nx-1));
        }
        else
        {
            for(i = 0; i <= nx-1; i++)
            {
                ptx(i) = 2*ap::randomreal()-1;
            }
        }
        
        //
        // Test:
        // * consistency of results from different queries
        // * points in query are IN the R-sphere (or at the boundary),
        //   and points not in query are outside of the R-sphere (or at the boundary)
        // * distances are correct and are ordered
        //
        kx = kdtreequeryaknn(treex, ptx, k, true, eps);
        kxy = kdtreequeryaknn(treexy, ptx, k, true, eps);
        kt = kdtreequeryaknn(treext, ptx, k, true, eps);
        if( kx!=k||kxy!=k||kt!=k )
        {
            kdterrors = true;
            return;
        }
        kx = 0;
        kxy = 0;
        kt = 0;
        kdtreequeryresultsx(treex, qx, kx);
        kdtreequeryresultsxy(treexy, qxy, kxy);
        kdtreequeryresultstags(treext, qtags, kt);
        kdtreequeryresultsdistances(treext, qr, kr);
        if( kx!=k||kxy!=k||kt!=k||kr!=k )
        {
            kdterrors = true;
            return;
        }
        kdterrors = kdterrors||kdtresultsdifferent(xy, n, qx, qxy, qtags, k, nx, ny);
        for(i = 0; i <= n-1; i++)
        {
            tmpb(i) = true;
        }
        r = 0;
        for(i = 0; i <= k-1; i++)
        {
            tmpb(qtags(i)) = false;
            ap::vmove(&tmpx(0), 1, &ptx(0), 1, ap::vlen(0,nx-1));
            ap::vsub(&tmpx(0), 1, &qx(i, 0), 1, ap::vlen(0,nx-1));
            r = ap::maxreal(r, vnorm(tmpx, nx, normtype));
        }
        for(i = 0; i <= n-1; i++)
        {
            if( tmpb(i) )
            {
                ap::vmove(&tmpx(0), 1, &ptx(0), 1, ap::vlen(0,nx-1));
                ap::vsub(&tmpx(0), 1, &xy(i, 0), 1, ap::vlen(0,nx-1));
                kdterrors = kdterrors||ap::fp_less(vnorm(tmpx, nx, normtype),r*(1-errtol)/(1+eps));
            }
        }
        for(i = 0; i <= k-2; i++)
        {
            kdterrors = kdterrors||ap::fp_greater(qr(i),qr(i+1));
        }
        for(i = 0; i <= k-1; i++)
        {
            ap::vmove(&tmpx(0), 1, &ptx(0), 1, ap::vlen(0,nx-1));
            ap::vsub(&tmpx(0), 1, &xy(qtags(i), 0), 1, ap::vlen(0,nx-1));
            kdterrors = kdterrors||ap::fp_greater(fabs(vnorm(tmpx, nx, normtype)-qr(i)),errtol);
        }
    }
    
    //
    // test general R-NN queries  (with self-matches):
    // * compare results from different trees (must be equal) and
    //   check that correct (value,tag) pairs are returned
    // * test results from XT tree - let R be radius of query result.
    //   then all points not in result must be not closer than R.
    //
    for(q = 1; q <= qcount; q++)
    {
        
        //
        // Select R
        //
        if( ap::fp_greater(ap::randomreal(),0.3) )
        {
            r = ap::maxreal(ap::randomreal(), ap::machineepsilon);
        }
        else
        {
            r = ap::machineepsilon;
        }
        
        //
        // Select point (either one of the points, or random)
        //
        if( ap::fp_greater(ap::randomreal(),0.5) )
        {
            i = ap::randominteger(n);
            ap::vmove(&ptx(0), 1, &xy(i, 0), 1, ap::vlen(0,nx-1));
        }
        else
        {
            for(i = 0; i <= nx-1; i++)
            {
                ptx(i) = 2*ap::randomreal()-1;
            }
        }
        
        //
        // Test:
        // * consistency of results from different queries
        // * points in query are IN the R-sphere (or at the boundary),
        //   and points not in query are outside of the R-sphere (or at the boundary)
        // * distances are correct and are ordered
        //
        kx = kdtreequeryrnn(treex, ptx, r, true);
        kxy = kdtreequeryrnn(treexy, ptx, r, true);
        kt = kdtreequeryrnn(treext, ptx, r, true);
        if( kxy!=kx||kt!=kx )
        {
            kdterrors = true;
            return;
        }
        kx = 0;
        kxy = 0;
        kt = 0;
        kdtreequeryresultsx(treex, qx, kx);
        kdtreequeryresultsxy(treexy, qxy, kxy);
        kdtreequeryresultstags(treext, qtags, kt);
        kdtreequeryresultsdistances(treext, qr, kr);
        if( kxy!=kx||kt!=kx||kr!=kx )
        {
            kdterrors = true;
            return;
        }
        kdterrors = kdterrors||kdtresultsdifferent(xy, n, qx, qxy, qtags, kx, nx, ny);
        for(i = 0; i <= n-1; i++)
        {
            tmpb(i) = true;
        }
        for(i = 0; i <= kx-1; i++)
        {
            tmpb(qtags(i)) = false;
        }
        for(i = 0; i <= n-1; i++)
        {
            ap::vmove(&tmpx(0), 1, &ptx(0), 1, ap::vlen(0,nx-1));
            ap::vsub(&tmpx(0), 1, &xy(i, 0), 1, ap::vlen(0,nx-1));
            if( tmpb(i) )
            {
                kdterrors = kdterrors||ap::fp_less(vnorm(tmpx, nx, normtype),r*(1-errtol));
            }
            else
            {
                kdterrors = kdterrors||ap::fp_greater(vnorm(tmpx, nx, normtype),r*(1+errtol));
            }
        }
        for(i = 0; i <= kx-2; i++)
        {
            kdterrors = kdterrors||ap::fp_greater(qr(i),qr(i+1));
        }
    }
    
    //
    // Test self-matching:
    // * self-match - nearest neighbor of each point in XY is the point itself
    // * no self-match - nearest neighbor is NOT the point itself
    //
    if( n>1 )
    {
        
        //
        // test for N=1 have non-general form, but it is not really needed
        //
        for(task = 0; task <= 1; task++)
        {
            for(i = 0; i <= n-1; i++)
            {
                ap::vmove(&ptx(0), 1, &xy(i, 0), 1, ap::vlen(0,nx-1));
                kx = kdtreequeryknn(treex, ptx, 1, task==0);
                kdtreequeryresultsx(treex, qx, kx);
                if( kx!=1 )
                {
                    kdterrors = true;
                    return;
                }
                isequal = true;
                for(j = 0; j <= nx-1; j++)
                {
                    isequal = isequal&&ap::fp_eq(qx(0,j),ptx(j));
                }
                if( task==0 )
                {
                    kdterrors = kdterrors||!isequal;
                }
                else
                {
                    kdterrors = kdterrors||isequal;
                }
            }
        }
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testnearestneighborunit_test_silent()
{
    bool result;

    result = testnearestneighbor(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testnearestneighborunit_test()
{
    bool result;

    result = testnearestneighbor(false);
    return result;
}





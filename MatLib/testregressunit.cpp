
 
#include <stdio.h>
#include "testregressunit.h"

static void generaterandomtask(double xl,
     double xr,
     bool randomx,
     double ymin,
     double ymax,
     double smin,
     double smax,
     int n,
     ap::real_2d_array& xy,
     ap::real_1d_array& s);
static void generatetask(double a,
     double b,
     double xl,
     double xr,
     bool randomx,
     double smin,
     double smax,
     int n,
     ap::real_2d_array& xy,
     ap::real_1d_array& s);
static void filltaskwithy(double a,
     double b,
     int n,
     ap::real_2d_array& xy,
     ap::real_1d_array& s);
static double generatenormal(double mean, double sigma);
static void calculatemv(const ap::real_1d_array& x,
     int n,
     double& mean,
     double& means,
     double& stddev,
     double& stddevs);
static void unsetlr(linearmodel& lr);

bool testlinregression(bool silent)
{
    bool result;
    double sigmathreshold;
    int maxn;
    int maxm;
    int passcount;
    int estpasscount;
    double threshold;
    int n;
    int i;
    int j;
    int k;
    int tmpi;
    int pass;
    int epass;
    int m;
    int tasktype;
    int modeltype;
    int m1;
    int m2;
    int n1;
    int n2;
    int info;
    int info2;
    ap::real_2d_array xy;
    ap::real_2d_array xy2;
    ap::real_1d_array s;
    ap::real_1d_array s2;
    ap::real_1d_array w2;
    ap::real_1d_array x;
    ap::real_1d_array ta;
    ap::real_1d_array tb;
    ap::real_1d_array tc;
    ap::real_1d_array xy0;
    ap::real_1d_array tmpweights;
    linearmodel w;
    linearmodel wt;
    linearmodel wt2;
    ap::real_1d_array x1;
    ap::real_1d_array x2;
    ap::real_1d_array ra;
    ap::real_1d_array ra2;
    double y1;
    double y2;
    int rlen;
    bool allsame;
    double ea;
    double eb;
    double varatested;
    double varbtested;
    double a;
    double b;
    double vara;
    double varb;
    double a2;
    double b2;
    double covab;
    double corrab;
    double p;
    int qcnt;
    ap::real_1d_array qtbl;
    ap::real_1d_array qvals;
    ap::real_1d_array qsigma;
    lrreport ar;
    lrreport ar2;
    double f;
    double fp;
    double fm;
    double v;
    double vv;
    double cvrmserror;
    double cvavgerror;
    double cvavgrelerror;
    double rmserror;
    double avgerror;
    double avgrelerror;
    bool nondefect;
    double sinshift;
    double tasklevel;
    double noiselevel;
    double hstep;
    double sigma;
    double mean;
    double means;
    double stddev;
    double stddevs;
    bool slcerrors;
    bool slerrors;
    bool grcoverrors;
    bool gropterrors;
    bool gresterrors;
    bool grothererrors;
    bool grconverrors;
    bool waserrors;

    
    //
    // Primary settings
    //
    maxn = 40;
    maxm = 5;
    passcount = 3;
    estpasscount = 1000;
    sigmathreshold = 7;
    threshold = 1000000*ap::machineepsilon;
    slerrors = false;
    slcerrors = false;
    grcoverrors = false;
    gropterrors = false;
    gresterrors = false;
    grothererrors = false;
    grconverrors = false;
    waserrors = false;
    
    //
    // Quantiles table setup
    //
    qcnt = 5;
    qtbl.setbounds(0, qcnt-1);
    qvals.setbounds(0, qcnt-1);
    qsigma.setbounds(0, qcnt-1);
    qtbl(0) = 0.5;
    qtbl(1) = 0.25;
    qtbl(2) = 0.10;
    qtbl(3) = 0.05;
    qtbl(4) = 0.025;
    for(i = 0; i <= qcnt-1; i++)
    {
        qsigma(i) = sqrt(qtbl(i)*(1-qtbl(i))/estpasscount);
    }
    
    //
    // Other setup
    //
    ta.setbounds(0, estpasscount-1);
    tb.setbounds(0, estpasscount-1);
    
    //
    // Test straight line regression
    //
    for(n = 2; n <= maxn; n++)
    {
        
        //
        // Fail/pass test
        //
        generaterandomtask(double(-1), double(1), false, double(-1), double(1), double(1), double(2), n, xy, s);
        lrlines(xy, s, n, info, a, b, vara, varb, covab, corrab, p);
        slcerrors = slcerrors||info!=1;
        generaterandomtask(double(+1), double(1), false, double(-1), double(1), double(1), double(2), n, xy, s);
        lrlines(xy, s, n, info, a, b, vara, varb, covab, corrab, p);
        slcerrors = slcerrors||info!=-3;
        generaterandomtask(double(-1), double(1), false, double(-1), double(1), double(-1), double(-1), n, xy, s);
        lrlines(xy, s, n, info, a, b, vara, varb, covab, corrab, p);
        slcerrors = slcerrors||info!=-2;
        generaterandomtask(double(-1), double(1), false, double(-1), double(1), double(2), double(1), 2, xy, s);
        lrlines(xy, s, 1, info, a, b, vara, varb, covab, corrab, p);
        slcerrors = slcerrors||info!=-1;
        
        //
        // Multipass tests
        //
        for(pass = 1; pass <= passcount; pass++)
        {
            
            //
            // Test S variant against non-S variant
            //
            ea = 2*ap::randomreal()-1;
            eb = 2*ap::randomreal()-1;
            generatetask(ea, eb, -5*ap::randomreal(), +5*ap::randomreal(), ap::fp_greater(ap::randomreal(),0.5), double(1), double(1), n, xy, s);
            lrlines(xy, s, n, info, a, b, vara, varb, covab, corrab, p);
            lrline(xy, n, info2, a2, b2);
            if( info!=1||info2!=1 )
            {
                slcerrors = true;
            }
            else
            {
                slerrors = slerrors||ap::fp_greater(fabs(a-a2),threshold)||ap::fp_greater(fabs(b-b2),threshold);
            }
            
            //
            // Test for A/B
            //
            // Generate task with exact, non-perturbed y[i],
            // then make non-zero s[i]
            //
            ea = 2*ap::randomreal()-1;
            eb = 2*ap::randomreal()-1;
            generatetask(ea, eb, -5*ap::randomreal(), +5*ap::randomreal(), n>4, 0.0, 0.0, n, xy, s);
            for(i = 0; i <= n-1; i++)
            {
                s(i) = 1+ap::randomreal();
            }
            lrlines(xy, s, n, info, a, b, vara, varb, covab, corrab, p);
            if( info!=1 )
            {
                slcerrors = true;
            }
            else
            {
                slerrors = slerrors||ap::fp_greater(fabs(a-ea),0.001)||ap::fp_greater(fabs(b-eb),0.001);
            }
            
            //
            // Test for VarA, VarB, P (P is being tested only for N>2)
            //
            for(i = 0; i <= qcnt-1; i++)
            {
                qvals(i) = 0;
            }
            ea = 2*ap::randomreal()-1;
            eb = 2*ap::randomreal()-1;
            generatetask(ea, eb, -5*ap::randomreal(), +5*ap::randomreal(), n>4, 1.0, 2.0, n, xy, s);
            lrlines(xy, s, n, info, a, b, vara, varb, covab, corrab, p);
            if( info!=1 )
            {
                slcerrors = true;
                continue;
            }
            varatested = vara;
            varbtested = varb;
            for(epass = 0; epass <= estpasscount-1; epass++)
            {
                
                //
                // Generate
                //
                filltaskwithy(ea, eb, n, xy, s);
                lrlines(xy, s, n, info, a, b, vara, varb, covab, corrab, p);
                if( info!=1 )
                {
                    slcerrors = true;
                    continue;
                }
                
                //
                // A, B, P
                // (P is being tested for uniformity, additional p-tests are below)
                //
                ta(epass) = a;
                tb(epass) = b;
                for(i = 0; i <= qcnt-1; i++)
                {
                    if( ap::fp_less_eq(p,qtbl(i)) )
                    {
                        qvals(i) = qvals(i)+double(1)/double(estpasscount);
                    }
                }
            }
            calculatemv(ta, estpasscount, mean, means, stddev, stddevs);
            slerrors = slerrors||ap::fp_greater_eq(fabs(mean-ea)/means,sigmathreshold);
            slerrors = slerrors||ap::fp_greater_eq(fabs(stddev-sqrt(varatested))/stddevs,sigmathreshold);
            calculatemv(tb, estpasscount, mean, means, stddev, stddevs);
            slerrors = slerrors||ap::fp_greater_eq(fabs(mean-eb)/means,sigmathreshold);
            slerrors = slerrors||ap::fp_greater_eq(fabs(stddev-sqrt(varbtested))/stddevs,sigmathreshold);
            if( n>2 )
            {
                for(i = 0; i <= qcnt-1; i++)
                {
                    if( ap::fp_greater(fabs(qtbl(i)-qvals(i))/qsigma(i),sigmathreshold) )
                    {
                        slerrors = true;
                    }
                }
            }
            
            //
            // Additional tests for P: correlation with fit quality
            //
            if( n>2 )
            {
                generatetask(ea, eb, -5*ap::randomreal(), +5*ap::randomreal(), false, 0.0, 0.0, n, xy, s);
                for(i = 0; i <= n-1; i++)
                {
                    s(i) = 1+ap::randomreal();
                }
                lrlines(xy, s, n, info, a, b, vara, varb, covab, corrab, p);
                if( info!=1 )
                {
                    slcerrors = true;
                    continue;
                }
                slerrors = slerrors||ap::fp_less(p,0.999);
                generatetask(double(0), double(0), -5*ap::randomreal(), +5*ap::randomreal(), false, 1.0, 1.0, n, xy, s);
                for(i = 0; i <= n-1; i++)
                {
                    if( i%2==0 )
                    {
                        xy(i,1) = +5.0;
                    }
                    else
                    {
                        xy(i,1) = -5.0;
                    }
                }
                if( n%2!=0 )
                {
                    xy(n-1,1) = 0;
                }
                lrlines(xy, s, n, info, a, b, vara, varb, covab, corrab, p);
                if( info!=1 )
                {
                    slcerrors = true;
                    continue;
                }
                slerrors = slerrors||ap::fp_greater(p,0.001);
            }
        }
    }
    
    //
    // General regression tests:
    //
    
    //
    // Simple linear tests (small sample, optimum point, covariance)
    //
    for(n = 3; n <= maxn; n++)
    {
        s.setbounds(0, n-1);
        
        //
        // Linear tests:
        // a. random points, sigmas
        // b. no sigmas
        //
        xy.setbounds(0, n-1, 0, 1);
        for(i = 0; i <= n-1; i++)
        {
            xy(i,0) = 2*ap::randomreal()-1;
            xy(i,1) = 2*ap::randomreal()-1;
            s(i) = 1+ap::randomreal();
        }
        lrbuilds(xy, s, n, 1, info, wt, ar);
        if( info!=1 )
        {
            grconverrors = true;
            continue;
        }
        lrunpack(wt, tmpweights, tmpi);
        lrlines(xy, s, n, info2, a, b, vara, varb, covab, corrab, p);
        gropterrors = gropterrors||ap::fp_greater(fabs(a-tmpweights(1)),threshold);
        gropterrors = gropterrors||ap::fp_greater(fabs(b-tmpweights(0)),threshold);
        grcoverrors = grcoverrors||ap::fp_greater(fabs(vara-ar.c(1,1)),threshold);
        grcoverrors = grcoverrors||ap::fp_greater(fabs(varb-ar.c(0,0)),threshold);
        grcoverrors = grcoverrors||ap::fp_greater(fabs(covab-ar.c(1,0)),threshold);
        grcoverrors = grcoverrors||ap::fp_greater(fabs(covab-ar.c(0,1)),threshold);
        lrbuild(xy, n, 1, info, wt, ar);
        if( info!=1 )
        {
            grconverrors = true;
            continue;
        }
        lrunpack(wt, tmpweights, tmpi);
        lrline(xy, n, info2, a, b);
        gropterrors = gropterrors||ap::fp_greater(fabs(a-tmpweights(1)),threshold);
        gropterrors = gropterrors||ap::fp_greater(fabs(b-tmpweights(0)),threshold);
    }
    
    //
    // S covariance versus S-less covariance.
    // Slightly skewed task, large sample size.
    // Will S-less subroutine estimate covariance matrix good enough?
    //
    n = 1000+ap::randominteger(3000);
    sigma = 0.1+ap::randomreal()*1.9;
    xy.setbounds(0, n-1, 0, 1);
    s.setbounds(0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        xy(i,0) = 1.5*ap::randomreal()-0.5;
        xy(i,1) = 1.2*xy(i,0)-0.3+generatenormal(double(0), sigma);
        s(i) = sigma;
    }
    lrbuild(xy, n, 1, info, wt, ar);
    lrlines(xy, s, n, info2, a, b, vara, varb, covab, corrab, p);
    if( info!=1||info2!=1 )
    {
        grconverrors = true;
    }
    else
    {
        grcoverrors = grcoverrors||ap::fp_greater(fabs(log(ar.c(0,0)/varb)),log(1.2));
        grcoverrors = grcoverrors||ap::fp_greater(fabs(log(ar.c(1,1)/vara)),log(1.2));
        grcoverrors = grcoverrors||ap::fp_greater(fabs(log(ar.c(0,1)/covab)),log(1.2));
        grcoverrors = grcoverrors||ap::fp_greater(fabs(log(ar.c(1,0)/covab)),log(1.2));
    }
    
    //
    // General tests:
    // * basis functions - up to cubic
    // * task types:
    // * data set is noisy sine half-period with random shift
    // * tests:
    //   unpacking/packing
    //   optimality
    //   error estimates
    // * tasks:
    //   0 = noised sine
    //   1 = degenerate task with 1-of-n encoded categorical variables
    //   2 = random task with large variation (for 1-type models)
    //   3 = random task with small variation (for 1-type models)
    //
    //   Additional tasks TODO
    //   specially designed task with defective vectors which leads to
    //   the failure of the fast CV formula.
    //
    //
    for(modeltype = 0; modeltype <= 1; modeltype++)
    {
        for(tasktype = 0; tasktype <= 3; tasktype++)
        {
            if( tasktype==0 )
            {
                m1 = 1;
                m2 = 3;
            }
            if( tasktype==1 )
            {
                m1 = 9;
                m2 = 9;
            }
            if( tasktype==2||tasktype==3 )
            {
                m1 = 9;
                m2 = 9;
            }
            for(m = m1; m <= m2; m++)
            {
                if( tasktype==0 )
                {
                    n1 = m+3;
                    n2 = m+20;
                }
                if( tasktype==1 )
                {
                    n1 = 70+ap::randominteger(70);
                    n2 = n1;
                }
                if( tasktype==2||tasktype==3 )
                {
                    n1 = 100;
                    n2 = n1;
                }
                for(n = n1; n <= n2; n++)
                {
                    xy.setbounds(0, n-1, 0, m);
                    xy0.setbounds(0, n-1);
                    s.setbounds(0, n-1);
                    hstep = 0.001;
                    noiselevel = 0.2;
                    
                    //
                    // Prepare task
                    //
                    if( tasktype==0 )
                    {
                        for(i = 0; i <= n-1; i++)
                        {
                            xy(i,0) = 2*ap::randomreal()-1;
                        }
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 1; j <= m-1; j++)
                            {
                                xy(i,j) = xy(i,0)*xy(i,j-1);
                            }
                        }
                        sinshift = ap::randomreal()*ap::pi();
                        for(i = 0; i <= n-1; i++)
                        {
                            xy0(i) = sin(sinshift+ap::pi()*0.5*(xy(i,0)+1));
                            xy(i,m) = xy0(i)+noiselevel*generatenormal(double(0), double(1));
                        }
                    }
                    if( tasktype==1 )
                    {
                        ap::ap_error::make_assertion(m==9, "");
                        ta.setbounds(0, 8);
                        ta(0) = 1;
                        ta(1) = 2;
                        ta(2) = 3;
                        ta(3) = 0.25;
                        ta(4) = 0.5;
                        ta(5) = 0.75;
                        ta(6) = 0.06;
                        ta(7) = 0.12;
                        ta(8) = 0.18;
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= m-1; j++)
                            {
                                xy(i,j) = 0;
                            }
                            xy(i,0+i%3) = 1;
                            xy(i,3+i/3%3) = 1;
                            xy(i,6+i/9%3) = 1;
                            v = ap::vdotproduct(&xy(i, 0), 1, &ta(0), 1, ap::vlen(0,8));
                            xy0(i) = v;
                            xy(i,m) = v+noiselevel*generatenormal(double(0), double(1));
                        }
                    }
                    if( tasktype==2||tasktype==3 )
                    {
                        ap::ap_error::make_assertion(m==9, "");
                        ta.setbounds(0, 8);
                        ta(0) = 1;
                        ta(1) = -2;
                        ta(2) = 3;
                        ta(3) = 0.25;
                        ta(4) = -0.5;
                        ta(5) = 0.75;
                        ta(6) = -0.06;
                        ta(7) = 0.12;
                        ta(8) = -0.18;
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= m-1; j++)
                            {
                                if( tasktype==2 )
                                {
                                    xy(i,j) = 1+generatenormal(double(0), double(3));
                                }
                                else
                                {
                                    xy(i,j) = 1+generatenormal(double(0), 0.05);
                                }
                            }
                            v = ap::vdotproduct(&xy(i, 0), 1, &ta(0), 1, ap::vlen(0,8));
                            xy0(i) = v;
                            xy(i,m) = v+noiselevel*generatenormal(double(0), double(1));
                        }
                    }
                    for(i = 0; i <= n-1; i++)
                    {
                        s(i) = 1+ap::randomreal();
                    }
                    
                    //
                    // Solve (using S-variant, non-S-variant is not tested)
                    //
                    if( modeltype==0 )
                    {
                        lrbuilds(xy, s, n, m, info, wt, ar);
                    }
                    else
                    {
                        lrbuildzs(xy, s, n, m, info, wt, ar);
                    }
                    if( info!=1 )
                    {
                        grconverrors = true;
                        continue;
                    }
                    lrunpack(wt, tmpweights, tmpi);
                    
                    //
                    // LRProcess test
                    //
                    x.setbounds(0, m-1);
                    v = tmpweights(m);
                    for(i = 0; i <= m-1; i++)
                    {
                        x(i) = 2*ap::randomreal()-1;
                        v = v+tmpweights(i)*x(i);
                    }
                    grothererrors = grothererrors||ap::fp_greater(fabs(v-lrprocess(wt, x))/ap::maxreal(fabs(v), double(1)),threshold);
                    
                    //
                    // LRPack test
                    //
                    lrpack(tmpweights, m, wt2);
                    x.setbounds(0, m-1);
                    for(i = 0; i <= m-1; i++)
                    {
                        x(i) = 2*ap::randomreal()-1;
                    }
                    v = lrprocess(wt, x);
                    grothererrors = grothererrors||ap::fp_greater(fabs(v-lrprocess(wt2, x))/fabs(v),threshold);
                    
                    //
                    // Optimality test
                    //
                    for(k = 0; k <= m; k++)
                    {
                        if( modeltype==1&&k==m )
                        {
                            
                            //
                            // 0-type models (with non-zero constant term)
                            // are tested for optimality of all coefficients.
                            //
                            // 1-type models (with zero constant term)
                            // are tested for optimality of non-constant terms only.
                            //
                            continue;
                        }
                        f = 0;
                        fp = 0;
                        fm = 0;
                        for(i = 0; i <= n-1; i++)
                        {
                            v = tmpweights(m);
                            for(j = 0; j <= m-1; j++)
                            {
                                v = v+xy(i,j)*tmpweights(j);
                            }
                            f = f+ap::sqr((v-xy(i,m))/s(i));
                            if( k<m )
                            {
                                vv = xy(i,k);
                            }
                            else
                            {
                                vv = 1;
                            }
                            fp = fp+ap::sqr((v+vv*hstep-xy(i,m))/s(i));
                            fm = fm+ap::sqr((v-vv*hstep-xy(i,m))/s(i));
                        }
                        gropterrors = gropterrors||ap::fp_greater(f,fp)||ap::fp_greater(f,fm);
                    }
                    
                    //
                    // Covariance matrix test:
                    // generate random vector, project coefficients on it,
                    // compare variance of projection with estimate provided
                    // by cov.matrix
                    //
                    ta.setbounds(0, estpasscount-1);
                    tb.setbounds(0, m);
                    tc.setbounds(0, m);
                    xy2.setbounds(0, n-1, 0, m);
                    for(i = 0; i <= m; i++)
                    {
                        tb(i) = generatenormal(double(0), double(1));
                    }
                    for(epass = 0; epass <= estpasscount-1; epass++)
                    {
                        for(i = 0; i <= n-1; i++)
                        {
                            ap::vmove(&xy2(i, 0), 1, &xy(i, 0), 1, ap::vlen(0,m-1));
                            xy2(i,m) = xy0(i)+s(i)*generatenormal(double(0), double(1));
                        }
                        if( modeltype==0 )
                        {
                            lrbuilds(xy2, s, n, m, info, wt, ar2);
                        }
                        else
                        {
                            lrbuildzs(xy2, s, n, m, info, wt, ar2);
                        }
                        if( info!=1 )
                        {
                            ta(epass) = 0;
                            grconverrors = true;
                            return result;
                        }
                        lrunpack(wt, w2, tmpi);
                        v = ap::vdotproduct(&tb(0), 1, &w2(0), 1, ap::vlen(0,m));
                        ta(epass) = v;
                    }
                    calculatemv(ta, estpasscount, mean, means, stddev, stddevs);
                    for(i = 0; i <= m; i++)
                    {
                        v = ap::vdotproduct(&tb(0), 1, &ar.c(0, i), ar.c.getstride(), ap::vlen(0,m));
                        tc(i) = v;
                    }
                    v = ap::vdotproduct(&tc(0), 1, &tb(0), 1, ap::vlen(0,m));
                    grcoverrors = grcoverrors||ap::fp_greater_eq(fabs((sqrt(v)-stddev)/stddevs),sigmathreshold);
                    
                    //
                    // Test for the fast CV error:
                    // calculate CV error by definition (leaving out N
                    // points and recalculating solution).
                    //
                    // Test for the training set error
                    //
                    cvrmserror = 0;
                    cvavgerror = 0;
                    cvavgrelerror = 0;
                    rmserror = 0;
                    avgerror = 0;
                    avgrelerror = 0;
                    xy2.setbounds(0, n-2, 0, m);
                    s2.setbounds(0, n-2);
                    for(i = 0; i <= n-2; i++)
                    {
                        ap::vmove(&xy2(i, 0), 1, &xy(i+1, 0), 1, ap::vlen(0,m));
                        s2(i) = s(i+1);
                    }
                    for(i = 0; i <= n-1; i++)
                    {
                        
                        //
                        // Trn
                        //
                        v = ap::vdotproduct(&xy(i, 0), 1, &tmpweights(0), 1, ap::vlen(0,m-1));
                        v = v+tmpweights(m);
                        rmserror = rmserror+ap::sqr(v-xy(i,m));
                        avgerror = avgerror+fabs(v-xy(i,m));
                        avgrelerror = avgrelerror+fabs((v-xy(i,m))/xy(i,m));
                        
                        //
                        // CV: non-defect vectors only
                        //
                        nondefect = true;
                        for(k = 0; k <= ar.ncvdefects-1; k++)
                        {
                            if( ar.cvdefects(k)==i )
                            {
                                nondefect = false;
                            }
                        }
                        if( nondefect )
                        {
                            if( modeltype==0 )
                            {
                                lrbuilds(xy2, s2, n-1, m, info2, wt, ar2);
                            }
                            else
                            {
                                lrbuildzs(xy2, s2, n-1, m, info2, wt, ar2);
                            }
                            if( info2!=1 )
                            {
                                grconverrors = true;
                                continue;
                            }
                            lrunpack(wt, w2, tmpi);
                            v = ap::vdotproduct(&xy(i, 0), 1, &w2(0), 1, ap::vlen(0,m-1));
                            v = v+w2(m);
                            cvrmserror = cvrmserror+ap::sqr(v-xy(i,m));
                            cvavgerror = cvavgerror+fabs(v-xy(i,m));
                            cvavgrelerror = cvavgrelerror+fabs((v-xy(i,m))/xy(i,m));
                        }
                        
                        //
                        // Next set
                        //
                        if( i!=n-1 )
                        {
                            ap::vmove(&xy2(i, 0), 1, &xy(i, 0), 1, ap::vlen(0,m));
                            s2(i) = s(i);
                        }
                    }
                    cvrmserror = sqrt(cvrmserror/(n-ar.ncvdefects));
                    cvavgerror = cvavgerror/(n-ar.ncvdefects);
                    cvavgrelerror = cvavgrelerror/(n-ar.ncvdefects);
                    rmserror = sqrt(rmserror/n);
                    avgerror = avgerror/n;
                    avgrelerror = avgrelerror/n;
                    gresterrors = gresterrors||ap::fp_greater(fabs(log(ar.cvrmserror/cvrmserror)),log(1+1.0E-5));
                    gresterrors = gresterrors||ap::fp_greater(fabs(log(ar.cvavgerror/cvavgerror)),log(1+1.0E-5));
                    gresterrors = gresterrors||ap::fp_greater(fabs(log(ar.cvavgrelerror/cvavgrelerror)),log(1+1.0E-5));
                    gresterrors = gresterrors||ap::fp_greater(fabs(log(ar.rmserror/rmserror)),log(1+1.0E-5));
                    gresterrors = gresterrors||ap::fp_greater(fabs(log(ar.avgerror/avgerror)),log(1+1.0E-5));
                    gresterrors = gresterrors||ap::fp_greater(fabs(log(ar.avgrelerror/avgrelerror)),log(1+1.0E-5));
                }
            }
        }
    }
    
    //
    // Additional subroutines
    //
    for(pass = 1; pass <= 50; pass++)
    {
        n = 2;
        do
        {
            noiselevel = ap::randomreal()+0.1;
            tasklevel = 2*ap::randomreal()-1;
        }
        while(ap::fp_less_eq(fabs(noiselevel-tasklevel),0.05));
        xy.setbounds(0, 3*n-1, 0, 1);
        for(i = 0; i <= n-1; i++)
        {
            xy(3*i+0,0) = i;
            xy(3*i+1,0) = i;
            xy(3*i+2,0) = i;
            xy(3*i+0,1) = tasklevel-noiselevel;
            xy(3*i+1,1) = tasklevel;
            xy(3*i+2,1) = tasklevel+noiselevel;
        }
        lrbuild(xy, 3*n, 1, info, wt, ar);
        if( info==1 )
        {
            lrunpack(wt, tmpweights, tmpi);
            v = lrrmserror(wt, xy, 3*n);
            grothererrors = grothererrors||ap::fp_greater(fabs(v-noiselevel*sqrt(double(2)/double(3))),threshold);
            v = lravgerror(wt, xy, 3*n);
            grothererrors = grothererrors||ap::fp_greater(fabs(v-noiselevel*(double(2)/double(3))),threshold);
            v = lravgrelerror(wt, xy, 3*n);
            vv = (fabs(noiselevel/(tasklevel-noiselevel))+fabs(noiselevel/(tasklevel+noiselevel)))/3;
            grothererrors = grothererrors||ap::fp_greater(fabs(v-vv),threshold*vv);
        }
        else
        {
            grothererrors = true;
        }
        for(i = 0; i <= n-1; i++)
        {
            xy(3*i+0,0) = i;
            xy(3*i+1,0) = i;
            xy(3*i+2,0) = i;
            xy(3*i+0,1) = -noiselevel;
            xy(3*i+1,1) = 0;
            xy(3*i+2,1) = +noiselevel;
        }
        lrbuild(xy, 3*n, 1, info, wt, ar);
        if( info==1 )
        {
            lrunpack(wt, tmpweights, tmpi);
            v = lravgrelerror(wt, xy, 3*n);
            grothererrors = grothererrors||ap::fp_greater(fabs(v-1),threshold);
        }
        else
        {
            grothererrors = true;
        }
    }
    for(pass = 1; pass <= 10; pass++)
    {
        m = 1+ap::randominteger(5);
        n = 10+ap::randominteger(10);
        xy.setbounds(0, n-1, 0, m);
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m; j++)
            {
                xy(i,j) = 2*ap::randomreal()-1;
            }
        }
        lrbuild(xy, n, m, info, w, ar);
        if( info<0 )
        {
            grothererrors = true;
            break;
        }
        x1.setbounds(0, m-1);
        x2.setbounds(0, m-1);
        
        //
        // Same inputs on original leads to same outputs
        // on copy created using LRCopy
        //
        unsetlr(wt);
        lrcopy(w, wt);
        for(i = 0; i <= m-1; i++)
        {
            x1(i) = 2*ap::randomreal()-1;
            x2(i) = x1(i);
        }
        y1 = lrprocess(w, x1);
        y2 = lrprocess(wt, x2);
        allsame = ap::fp_eq(y1,y2);
        grothererrors = grothererrors||!allsame;
        
        //
        // Same inputs on original leads to same outputs
        // on copy created using LRSerialize
        //
        unsetlr(wt);
        ra.setbounds(0, 0);
        ra(0) = 0;
        rlen = 0;
        lrserialize(w, ra, rlen);
        ra2.setbounds(0, rlen-1);
        for(i = 0; i <= rlen-1; i++)
        {
            ra2(i) = ra(i);
        }
        lrunserialize(ra2, wt);
        for(i = 0; i <= m-1; i++)
        {
            x1(i) = 2*ap::randomreal()-1;
            x2(i) = x1(i);
        }
        y1 = lrprocess(w, x1);
        y2 = lrprocess(wt, x2);
        allsame = ap::fp_eq(y1,y2);
        grothererrors = grothererrors||!allsame;
    }
    
    //
    // TODO: Degenerate tests (when design matrix and right part are zero)
    //
    
    //
    // Final report
    //
    waserrors = slerrors||slcerrors||gropterrors||grcoverrors||gresterrors||grothererrors||grconverrors;
    if( !silent )
    {
        printf("REGRESSION TEST\n");
        printf("STRAIGHT LINE REGRESSION:                ");
        if( !slerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("STRAIGHT LINE REGRESSION CONVERGENCE:    ");
        if( !slcerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("GENERAL LINEAR REGRESSION:               ");
        if( !(gropterrors||grcoverrors||gresterrors||grothererrors||grconverrors) )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* OPTIMALITY:                            ");
        if( !gropterrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* COV. MATRIX:                           ");
        if( !grcoverrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* ERROR ESTIMATES:                       ");
        if( !gresterrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* CONVERGENCE:                           ");
        if( !grconverrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* OTHER SUBROUTINES:                     ");
        if( !grothererrors )
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
Task generation. Meaningless task, just random numbers.
*************************************************************************/
static void generaterandomtask(double xl,
     double xr,
     bool randomx,
     double ymin,
     double ymax,
     double smin,
     double smax,
     int n,
     ap::real_2d_array& xy,
     ap::real_1d_array& s)
{
    int i;

    xy.setbounds(0, n-1, 0, 1);
    s.setbounds(0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        if( randomx )
        {
            xy(i,0) = xl+(xr-xl)*ap::randomreal();
        }
        else
        {
            xy(i,0) = xl+(xr-xl)*i/(n-1);
        }
        xy(i,1) = ymin+(ymax-ymin)*ap::randomreal();
        s(i) = smin+(smax-smin)*ap::randomreal();
    }
}


/*************************************************************************
Task generation.
*************************************************************************/
static void generatetask(double a,
     double b,
     double xl,
     double xr,
     bool randomx,
     double smin,
     double smax,
     int n,
     ap::real_2d_array& xy,
     ap::real_1d_array& s)
{
    int i;

    xy.setbounds(0, n-1, 0, 1);
    s.setbounds(0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        if( randomx )
        {
            xy(i,0) = xl+(xr-xl)*ap::randomreal();
        }
        else
        {
            xy(i,0) = xl+(xr-xl)*i/(n-1);
        }
        s(i) = smin+(smax-smin)*ap::randomreal();
        xy(i,1) = a+b*xy(i,0)+generatenormal(double(0), s(i));
    }
}


/*************************************************************************
Task generation.
y[i] are filled based on A, B, X[I], S[I]
*************************************************************************/
static void filltaskwithy(double a,
     double b,
     int n,
     ap::real_2d_array& xy,
     ap::real_1d_array& s)
{
    int i;

    for(i = 0; i <= n-1; i++)
    {
        xy(i,1) = a+b*xy(i,0)+generatenormal(double(0), s(i));
    }
}


/*************************************************************************
Normal random numbers
*************************************************************************/
static double generatenormal(double mean, double sigma)
{
    double result;
    double u;
    double v;
    double s;
    double sum;

    result = mean;
    while(true)
    {
        u = (2*ap::randominteger(2)-1)*ap::randomreal();
        v = (2*ap::randominteger(2)-1)*ap::randomreal();
        sum = u*u+v*v;
        if( ap::fp_less(sum,1)&&ap::fp_greater(sum,0) )
        {
            sum = sqrt(-2*log(sum)/sum);
            result = sigma*u*sum+mean;
            return result;
        }
    }
    return result;
}


/*************************************************************************
Moments estimates and their errors
*************************************************************************/
static void calculatemv(const ap::real_1d_array& x,
     int n,
     double& mean,
     double& means,
     double& stddev,
     double& stddevs)
{
    int i;
    double v;
    double v1;
    double v2;
    double variance;

    mean = 0;
    means = 1;
    stddev = 0;
    stddevs = 1;
    variance = 0;
    if( n<=1 )
    {
        return;
    }
    
    //
    // Mean
    //
    for(i = 0; i <= n-1; i++)
    {
        mean = mean+x(i);
    }
    mean = mean/n;
    
    //
    // Variance (using corrected two-pass algorithm)
    //
    if( n!=1 )
    {
        v1 = 0;
        for(i = 0; i <= n-1; i++)
        {
            v1 = v1+ap::sqr(x(i)-mean);
        }
        v2 = 0;
        for(i = 0; i <= n-1; i++)
        {
            v2 = v2+(x(i)-mean);
        }
        v2 = ap::sqr(v2)/n;
        variance = (v1-v2)/(n-1);
        if( ap::fp_less(variance,0) )
        {
            variance = 0;
        }
        stddev = sqrt(variance);
    }
    
    //
    // Errors
    //
    means = stddev/sqrt(double(n));
    stddevs = stddev*sqrt(double(2))/sqrt(double(n-1));
}


/*************************************************************************
Unsets LR
*************************************************************************/
static void unsetlr(linearmodel& lr)
{
    ap::real_2d_array xy;
    int info;
    lrreport rep;
    int i;

    xy.setbounds(0, 5, 0, 1);
    for(i = 0; i <= 5; i++)
    {
        xy(i,0) = 0;
        xy(i,1) = 0;
    }
    lrbuild(xy, 6, 1, info, lr, rep);
    ap::ap_error::make_assertion(info>0, "");
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testregressunit_test_silent()
{
    bool result;

    result = testlinregression(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testregressunit_test()
{
    bool result;

    result = testlinregression(false);
    return result;
}





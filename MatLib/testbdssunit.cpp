
 
#include <stdio.h>
#include "testbdssunit.h"

static void unset2d(ap::complex_2d_array& a);
static void unset1d(ap::real_1d_array& a);
static void unset1di(ap::integer_1d_array& a);
static void testsortresults(const ap::real_1d_array& asorted,
     const ap::integer_1d_array& p1,
     const ap::integer_1d_array& p2,
     const ap::real_1d_array& aoriginal,
     int n,
     bool& waserrors);

/*************************************************************************
Testing BDSS operations
*************************************************************************/
bool testbdss(bool silent)
{
    bool result;
    int n;
    int i;
    int j;
    int pass;
    int passcount;
    int maxn;
    int maxnq;
    ap::real_1d_array a;
    ap::real_1d_array a0;
    ap::real_1d_array at;
    ap::real_2d_array p;
    ap::real_1d_array thresholds;
    int ni;
    ap::integer_1d_array c;
    ap::integer_1d_array p1;
    ap::integer_1d_array p2;
    ap::integer_1d_array ties;
    ap::integer_1d_array pt1;
    ap::integer_1d_array pt2;
    int tiecount;
    int c1;
    int c0;
    int nc;
    ap::real_1d_array tmp;
    double pal;
    double pbl;
    double par;
    double pbr;
    double cve;
    double cvr;
    int info;
    double threshold;
    ap::integer_1d_array tiebuf;
    ap::integer_1d_array cntbuf;
    double rms;
    double cvrms;
    bool waserrors;
    bool tieserrors;
    bool split2errors;
    bool optimalsplitkerrors;
    bool splitkerrors;

    waserrors = false;
    tieserrors = false;
    split2errors = false;
    splitkerrors = false;
    optimalsplitkerrors = false;
    maxn = 100;
    maxnq = 49;
    passcount = 10;
    
    //
    // Test ties
    //
    for(n = 1; n <= maxn; n++)
    {
        for(pass = 1; pass <= passcount; pass++)
        {
            
            //
            // untied data, test DSTie
            //
            unset1di(p1);
            unset1di(p2);
            unset1di(pt1);
            unset1di(pt2);
            a.setbounds(0, n-1);
            a0.setbounds(0, n-1);
            at.setbounds(0, n-1);
            tmp.setbounds(0, n-1);
            a(0) = 2*ap::randomreal()-1;
            tmp(0) = ap::randomreal();
            for(i = 1; i <= n-1; i++)
            {
                
                //
                // A is randomly permuted
                //
                a(i) = a(i-1)+0.1*ap::randomreal()+0.1;
                tmp(i) = ap::randomreal();
            }
            tagsortfastr(tmp, a, n);
            for(i = 0; i <= n-1; i++)
            {
                a0(i) = a(i);
                at(i) = a(i);
            }
            dstie(a0, n, ties, tiecount, p1, p2);
            tagsort(at, n, pt1, pt2);
            for(i = 0; i <= n-1; i++)
            {
                tieserrors = tieserrors||p1(i)!=pt1(i);
                tieserrors = tieserrors||p2(i)!=pt2(i);
            }
            tieserrors = tieserrors||tiecount!=n;
            if( tiecount==n )
            {
                for(i = 0; i <= n; i++)
                {
                    tieserrors = tieserrors||ties(i)!=i;
                }
            }
            
            //
            // tied data, test DSTie
            //
            unset1di(p1);
            unset1di(p2);
            unset1di(pt1);
            unset1di(pt2);
            a.setbounds(0, n-1);
            a0.setbounds(0, n-1);
            at.setbounds(0, n-1);
            c1 = 0;
            c0 = 0;
            for(i = 0; i <= n-1; i++)
            {
                a(i) = ap::randominteger(2);
                if( ap::fp_eq(a(i),0) )
                {
                    c0 = c0+1;
                }
                else
                {
                    c1 = c1+1;
                }
                a0(i) = a(i);
                at(i) = a(i);
            }
            dstie(a0, n, ties, tiecount, p1, p2);
            tagsort(at, n, pt1, pt2);
            for(i = 0; i <= n-1; i++)
            {
                tieserrors = tieserrors||p1(i)!=pt1(i);
                tieserrors = tieserrors||p2(i)!=pt2(i);
            }
            if( c0==0||c1==0 )
            {
                tieserrors = tieserrors||tiecount!=1;
                if( tiecount==1 )
                {
                    tieserrors = tieserrors||ties(0)!=0;
                    tieserrors = tieserrors||ties(1)!=n;
                }
            }
            else
            {
                tieserrors = tieserrors||tiecount!=2;
                if( tiecount==2 )
                {
                    tieserrors = tieserrors||ties(0)!=0;
                    tieserrors = tieserrors||ties(1)!=c0;
                    tieserrors = tieserrors||ties(2)!=n;
                }
            }
        }
    }
    
    //
    // split-2
    //
    
    //
    // General tests for different N's
    //
    for(n = 1; n <= maxn; n++)
    {
        a.setbounds(0, n-1);
        c.setbounds(0, n-1);
        
        //
        // one-tie test
        //
        if( n%2==0 )
        {
            for(i = 0; i <= n-1; i++)
            {
                a(i) = n;
                c(i) = i%2;
            }
            dsoptimalsplit2(a, c, n, info, threshold, pal, pbl, par, pbr, cve);
            if( info!=-3 )
            {
                split2errors = true;
                continue;
            }
        }
        
        //
        // two-tie test
        //
        
        //
        // test #1
        //
        if( n>1 )
        {
            for(i = 0; i <= n-1; i++)
            {
                a(i) = i/((n+1)/2);
                c(i) = i/((n+1)/2);
            }
            dsoptimalsplit2(a, c, n, info, threshold, pal, pbl, par, pbr, cve);
            if( info!=1 )
            {
                split2errors = true;
                continue;
            }
            split2errors = split2errors||ap::fp_greater(fabs(threshold-0.5),100*ap::machineepsilon);
            split2errors = split2errors||ap::fp_greater(fabs(pal-1),100*ap::machineepsilon);
            split2errors = split2errors||ap::fp_greater(fabs(pbl-0),100*ap::machineepsilon);
            split2errors = split2errors||ap::fp_greater(fabs(par-0),100*ap::machineepsilon);
            split2errors = split2errors||ap::fp_greater(fabs(pbr-1),100*ap::machineepsilon);
        }
    }
    
    //
    // Special "CREDIT"-test (transparency coefficient)
    //
    n = 110;
    a.setbounds(0, n-1);
    c.setbounds(0, n-1);
    a(0) = 0.000;
    c(0) = 0;
    a(1) = 0.000;
    c(1) = 0;
    a(2) = 0.000;
    c(2) = 0;
    a(3) = 0.000;
    c(3) = 0;
    a(4) = 0.000;
    c(4) = 0;
    a(5) = 0.000;
    c(5) = 0;
    a(6) = 0.000;
    c(6) = 0;
    a(7) = 0.000;
    c(7) = 1;
    a(8) = 0.000;
    c(8) = 0;
    a(9) = 0.000;
    c(9) = 1;
    a(10) = 0.000;
    c(10) = 0;
    a(11) = 0.000;
    c(11) = 0;
    a(12) = 0.000;
    c(12) = 0;
    a(13) = 0.000;
    c(13) = 0;
    a(14) = 0.000;
    c(14) = 0;
    a(15) = 0.000;
    c(15) = 0;
    a(16) = 0.000;
    c(16) = 0;
    a(17) = 0.000;
    c(17) = 0;
    a(18) = 0.000;
    c(18) = 0;
    a(19) = 0.000;
    c(19) = 0;
    a(20) = 0.000;
    c(20) = 0;
    a(21) = 0.000;
    c(21) = 0;
    a(22) = 0.000;
    c(22) = 1;
    a(23) = 0.000;
    c(23) = 0;
    a(24) = 0.000;
    c(24) = 0;
    a(25) = 0.000;
    c(25) = 0;
    a(26) = 0.000;
    c(26) = 0;
    a(27) = 0.000;
    c(27) = 1;
    a(28) = 0.000;
    c(28) = 0;
    a(29) = 0.000;
    c(29) = 1;
    a(30) = 0.000;
    c(30) = 0;
    a(31) = 0.000;
    c(31) = 1;
    a(32) = 0.000;
    c(32) = 0;
    a(33) = 0.000;
    c(33) = 1;
    a(34) = 0.000;
    c(34) = 0;
    a(35) = 0.030;
    c(35) = 0;
    a(36) = 0.030;
    c(36) = 0;
    a(37) = 0.050;
    c(37) = 0;
    a(38) = 0.070;
    c(38) = 1;
    a(39) = 0.110;
    c(39) = 0;
    a(40) = 0.110;
    c(40) = 1;
    a(41) = 0.120;
    c(41) = 0;
    a(42) = 0.130;
    c(42) = 0;
    a(43) = 0.140;
    c(43) = 0;
    a(44) = 0.140;
    c(44) = 0;
    a(45) = 0.140;
    c(45) = 0;
    a(46) = 0.150;
    c(46) = 0;
    a(47) = 0.150;
    c(47) = 0;
    a(48) = 0.170;
    c(48) = 0;
    a(49) = 0.190;
    c(49) = 1;
    a(50) = 0.200;
    c(50) = 0;
    a(51) = 0.200;
    c(51) = 0;
    a(52) = 0.250;
    c(52) = 0;
    a(53) = 0.250;
    c(53) = 0;
    a(54) = 0.260;
    c(54) = 0;
    a(55) = 0.270;
    c(55) = 0;
    a(56) = 0.280;
    c(56) = 0;
    a(57) = 0.310;
    c(57) = 0;
    a(58) = 0.310;
    c(58) = 0;
    a(59) = 0.330;
    c(59) = 0;
    a(60) = 0.330;
    c(60) = 0;
    a(61) = 0.340;
    c(61) = 0;
    a(62) = 0.340;
    c(62) = 0;
    a(63) = 0.370;
    c(63) = 0;
    a(64) = 0.380;
    c(64) = 1;
    a(65) = 0.380;
    c(65) = 0;
    a(66) = 0.410;
    c(66) = 0;
    a(67) = 0.460;
    c(67) = 0;
    a(68) = 0.520;
    c(68) = 0;
    a(69) = 0.530;
    c(69) = 0;
    a(70) = 0.540;
    c(70) = 0;
    a(71) = 0.560;
    c(71) = 0;
    a(72) = 0.560;
    c(72) = 0;
    a(73) = 0.570;
    c(73) = 0;
    a(74) = 0.600;
    c(74) = 0;
    a(75) = 0.600;
    c(75) = 0;
    a(76) = 0.620;
    c(76) = 0;
    a(77) = 0.650;
    c(77) = 0;
    a(78) = 0.660;
    c(78) = 0;
    a(79) = 0.680;
    c(79) = 0;
    a(80) = 0.700;
    c(80) = 0;
    a(81) = 0.750;
    c(81) = 0;
    a(82) = 0.770;
    c(82) = 0;
    a(83) = 0.770;
    c(83) = 0;
    a(84) = 0.770;
    c(84) = 0;
    a(85) = 0.790;
    c(85) = 0;
    a(86) = 0.810;
    c(86) = 0;
    a(87) = 0.840;
    c(87) = 0;
    a(88) = 0.860;
    c(88) = 0;
    a(89) = 0.870;
    c(89) = 0;
    a(90) = 0.890;
    c(90) = 0;
    a(91) = 0.900;
    c(91) = 1;
    a(92) = 0.900;
    c(92) = 0;
    a(93) = 0.910;
    c(93) = 0;
    a(94) = 0.940;
    c(94) = 0;
    a(95) = 0.950;
    c(95) = 0;
    a(96) = 0.952;
    c(96) = 0;
    a(97) = 0.970;
    c(97) = 0;
    a(98) = 0.970;
    c(98) = 0;
    a(99) = 0.980;
    c(99) = 0;
    a(100) = 1.000;
    c(100) = 0;
    a(101) = 1.000;
    c(101) = 0;
    a(102) = 1.000;
    c(102) = 0;
    a(103) = 1.000;
    c(103) = 0;
    a(104) = 1.000;
    c(104) = 0;
    a(105) = 1.020;
    c(105) = 0;
    a(106) = 1.090;
    c(106) = 0;
    a(107) = 1.130;
    c(107) = 0;
    a(108) = 1.840;
    c(108) = 0;
    a(109) = 2.470;
    c(109) = 0;
    dsoptimalsplit2(a, c, n, info, threshold, pal, pbl, par, pbr, cve);
    if( info!=1 )
    {
        split2errors = true;
    }
    else
    {
        split2errors = split2errors||ap::fp_greater(fabs(threshold-0.195),100*ap::machineepsilon);
        split2errors = split2errors||ap::fp_greater(fabs(pal-0.80),0.02);
        split2errors = split2errors||ap::fp_greater(fabs(pbl-0.20),0.02);
        split2errors = split2errors||ap::fp_greater(fabs(par-0.97),0.02);
        split2errors = split2errors||ap::fp_greater(fabs(pbr-0.03),0.02);
    }
    
    //
    // split-2 fast
    //
    
    //
    // General tests for different N's
    //
    for(n = 1; n <= maxn; n++)
    {
        a.setbounds(0, n-1);
        c.setbounds(0, n-1);
        tiebuf.setbounds(0, n);
        cntbuf.setbounds(0, 3);
        
        //
        // one-tie test
        //
        if( n%2==0 )
        {
            for(i = 0; i <= n-1; i++)
            {
                a(i) = n;
                c(i) = i%2;
            }
            dsoptimalsplit2fast(a, c, tiebuf, cntbuf, n, 2, 0.00, info, threshold, rms, cvrms);
            if( info!=-3 )
            {
                split2errors = true;
                continue;
            }
        }
        
        //
        // two-tie test
        //
        
        //
        // test #1
        //
        if( n>1 )
        {
            for(i = 0; i <= n-1; i++)
            {
                a(i) = i/((n+1)/2);
                c(i) = i/((n+1)/2);
            }
            dsoptimalsplit2fast(a, c, tiebuf, cntbuf, n, 2, 0.00, info, threshold, rms, cvrms);
            if( info!=1 )
            {
                split2errors = true;
                continue;
            }
            split2errors = split2errors||ap::fp_greater(fabs(threshold-0.5),100*ap::machineepsilon);
            split2errors = split2errors||ap::fp_greater(fabs(rms-0),100*ap::machineepsilon);
            if( n==2 )
            {
                split2errors = split2errors||ap::fp_greater(fabs(cvrms-0.5),100*ap::machineepsilon);
            }
            else
            {
                if( n==3 )
                {
                    split2errors = split2errors||ap::fp_greater(fabs(cvrms-sqrt((2*0+2*0+2*0.25)/6)),100*ap::machineepsilon);
                }
                else
                {
                    split2errors = split2errors||ap::fp_greater(fabs(cvrms),100*ap::machineepsilon);
                }
            }
        }
    }
    
    //
    // special tests
    //
    n = 10;
    a.setbounds(0, n-1);
    c.setbounds(0, n-1);
    tiebuf.setbounds(0, n);
    cntbuf.setbounds(0, 2*3-1);
    for(i = 0; i <= n-1; i++)
    {
        a(i) = i;
        if( i<=n-3 )
        {
            c(i) = 0;
        }
        else
        {
            c(i) = i-(n-3);
        }
    }
    dsoptimalsplit2fast(a, c, tiebuf, cntbuf, n, 3, 0.00, info, threshold, rms, cvrms);
    if( info!=1 )
    {
        split2errors = true;
    }
    else
    {
        split2errors = split2errors||ap::fp_greater(fabs(threshold-(n-2.5)),100*ap::machineepsilon);
        split2errors = split2errors||ap::fp_greater(fabs(rms-sqrt((0.25+0.25+0.25+0.25)/(3*n))),100*ap::machineepsilon);
        split2errors = split2errors||ap::fp_greater(fabs(cvrms-sqrt(double(1+1+1+1)/double(3*n))),100*ap::machineepsilon);
    }
    
    //
    // Optimal split-K
    //
    
    //
    // General tests for different N's
    //
    for(n = 1; n <= maxnq; n++)
    {
        a.setbounds(0, n-1);
        c.setbounds(0, n-1);
        
        //
        // one-tie test
        //
        if( n%2==0 )
        {
            for(i = 0; i <= n-1; i++)
            {
                a(i) = pass;
                c(i) = i%2;
            }
            dsoptimalsplitk(a, c, n, 2, 2+ap::randominteger(5), info, thresholds, ni, cve);
            if( info!=-3 )
            {
                optimalsplitkerrors = true;
                continue;
            }
        }
        
        //
        // two-tie test
        //
        
        //
        // test #1
        //
        if( n>1 )
        {
            c0 = 0;
            c1 = 0;
            for(i = 0; i <= n-1; i++)
            {
                a(i) = i/((n+1)/2);
                c(i) = i/((n+1)/2);
                if( c(i)==0 )
                {
                    c0 = c0+1;
                }
                if( c(i)==1 )
                {
                    c1 = c1+1;
                }
            }
            dsoptimalsplitk(a, c, n, 2, 2+ap::randominteger(5), info, thresholds, ni, cve);
            if( info!=1 )
            {
                optimalsplitkerrors = true;
                continue;
            }
            optimalsplitkerrors = optimalsplitkerrors||ni!=2;
            optimalsplitkerrors = optimalsplitkerrors||ap::fp_greater(fabs(thresholds(0)-0.5),100*ap::machineepsilon);
            optimalsplitkerrors = optimalsplitkerrors||ap::fp_greater(fabs(cve-(-c0*log(double(c0)/double(c0+1))-c1*log(double(c1)/double(c1+1)))),100*ap::machineepsilon);
        }
        
        //
        // test #2
        //
        if( n>2 )
        {
            c0 = 1+ap::randominteger(n-1);
            c1 = n-c0;
            for(i = 0; i <= n-1; i++)
            {
                if( i<c0 )
                {
                    a(i) = 0;
                    c(i) = 0;
                }
                else
                {
                    a(i) = 1;
                    c(i) = 1;
                }
            }
            dsoptimalsplitk(a, c, n, 2, 2+ap::randominteger(5), info, thresholds, ni, cve);
            if( info!=1 )
            {
                optimalsplitkerrors = true;
                continue;
            }
            optimalsplitkerrors = optimalsplitkerrors||ni!=2;
            optimalsplitkerrors = optimalsplitkerrors||ap::fp_greater(fabs(thresholds(0)-0.5),100*ap::machineepsilon);
            optimalsplitkerrors = optimalsplitkerrors||ap::fp_greater(fabs(cve-(-c0*log(double(c0)/double(c0+1))-c1*log(double(c1)/double(c1+1)))),100*ap::machineepsilon);
        }
        
        //
        // multi-tie test
        //
        if( n>=16 )
        {
            
            //
            // Multi-tie test.
            //
            // First NC-1 ties have C0 entries, remaining NC-th tie
            // have C1 entries.
            //
            nc = ap::round(sqrt(double(n)));
            c0 = n/nc;
            c1 = n-c0*(nc-1);
            for(i = 0; i <= nc-2; i++)
            {
                for(j = c0*i; j <= c0*(i+1)-1; j++)
                {
                    a(j) = j;
                    c(j) = i;
                }
            }
            for(j = c0*(nc-1); j <= n-1; j++)
            {
                a(j) = j;
                c(j) = nc-1;
            }
            dsoptimalsplitk(a, c, n, nc, nc+ap::randominteger(nc), info, thresholds, ni, cve);
            if( info!=1 )
            {
                optimalsplitkerrors = true;
                continue;
            }
            optimalsplitkerrors = optimalsplitkerrors||ni!=nc;
            if( ni==nc )
            {
                for(i = 0; i <= nc-2; i++)
                {
                    optimalsplitkerrors = optimalsplitkerrors||ap::fp_greater(fabs(thresholds(i)-(c0*(i+1)-1+0.5)),100*ap::machineepsilon);
                }
                cvr = -((nc-1)*c0*log(double(c0)/double(c0+nc-1))+c1*log(double(c1)/double(c1+nc-1)));
                optimalsplitkerrors = optimalsplitkerrors||ap::fp_greater(fabs(cve-cvr),100*ap::machineepsilon);
            }
        }
    }
    
    //
    // Non-optimal split-K
    //
    
    //
    // General tests for different N's
    //
    for(n = 1; n <= maxnq; n++)
    {
        a.setbounds(0, n-1);
        c.setbounds(0, n-1);
        
        //
        // one-tie test
        //
        if( n%2==0 )
        {
            for(i = 0; i <= n-1; i++)
            {
                a(i) = pass;
                c(i) = i%2;
            }
            dssplitk(a, c, n, 2, 2+ap::randominteger(5), info, thresholds, ni, cve);
            if( info!=-3 )
            {
                splitkerrors = true;
                continue;
            }
        }
        
        //
        // two-tie test
        //
        
        //
        // test #1
        //
        if( n>1 )
        {
            c0 = 0;
            c1 = 0;
            for(i = 0; i <= n-1; i++)
            {
                a(i) = i/((n+1)/2);
                c(i) = i/((n+1)/2);
                if( c(i)==0 )
                {
                    c0 = c0+1;
                }
                if( c(i)==1 )
                {
                    c1 = c1+1;
                }
            }
            dssplitk(a, c, n, 2, 2+ap::randominteger(5), info, thresholds, ni, cve);
            if( info!=1 )
            {
                splitkerrors = true;
                continue;
            }
            splitkerrors = splitkerrors||ni!=2;
            if( ni==2 )
            {
                splitkerrors = splitkerrors||ap::fp_greater(fabs(thresholds(0)-0.5),100*ap::machineepsilon);
                splitkerrors = splitkerrors||ap::fp_greater(fabs(cve-(-c0*log(double(c0)/double(c0+1))-c1*log(double(c1)/double(c1+1)))),100*ap::machineepsilon);
            }
        }
        
        //
        // test #2
        //
        if( n>2 )
        {
            c0 = 1+ap::randominteger(n-1);
            c1 = n-c0;
            for(i = 0; i <= n-1; i++)
            {
                if( i<c0 )
                {
                    a(i) = 0;
                    c(i) = 0;
                }
                else
                {
                    a(i) = 1;
                    c(i) = 1;
                }
            }
            dssplitk(a, c, n, 2, 2+ap::randominteger(5), info, thresholds, ni, cve);
            if( info!=1 )
            {
                splitkerrors = true;
                continue;
            }
            splitkerrors = splitkerrors||ni!=2;
            if( ni==2 )
            {
                splitkerrors = splitkerrors||ap::fp_greater(fabs(thresholds(0)-0.5),100*ap::machineepsilon);
                splitkerrors = splitkerrors||ap::fp_greater(fabs(cve-(-c0*log(double(c0)/double(c0+1))-c1*log(double(c1)/double(c1+1)))),100*ap::machineepsilon);
            }
        }
        
        //
        // multi-tie test
        //
        for(c0 = 4; c0 <= n; c0++)
        {
            if( n%c0==0&&n/c0<=c0&&n/c0>1 )
            {
                nc = n/c0;
                for(i = 0; i <= nc-1; i++)
                {
                    for(j = c0*i; j <= c0*(i+1)-1; j++)
                    {
                        a(j) = j;
                        c(j) = i;
                    }
                }
                dssplitk(a, c, n, nc, nc+ap::randominteger(nc), info, thresholds, ni, cve);
                if( info!=1 )
                {
                    splitkerrors = true;
                    continue;
                }
                splitkerrors = splitkerrors||ni!=nc;
                if( ni==nc )
                {
                    for(i = 0; i <= nc-2; i++)
                    {
                        splitkerrors = splitkerrors||ap::fp_greater(fabs(thresholds(i)-(c0*(i+1)-1+0.5)),100*ap::machineepsilon);
                    }
                    cvr = -nc*c0*log(double(c0)/double(c0+nc-1));
                    splitkerrors = splitkerrors||ap::fp_greater(fabs(cve-cvr),100*ap::machineepsilon);
                }
            }
        }
    }
    
    //
    // report
    //
    waserrors = tieserrors||split2errors||optimalsplitkerrors||splitkerrors;
    if( !silent )
    {
        printf("TESTING BASIC DATASET SUBROUTINES\n");
        printf("TIES:                               ");
        if( !tieserrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("SPLIT-2:                            ");
        if( !split2errors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("OPTIMAL SPLIT-K:                    ");
        if( !optimalsplitkerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("SPLIT-K:                            ");
        if( !splitkerrors )
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
Unsets 1D array.
*************************************************************************/
static void unset1di(ap::integer_1d_array& a)
{

    a.setbounds(0, 0);
    a(0) = ap::randominteger(3)-1;
}


static void testsortresults(const ap::real_1d_array& asorted,
     const ap::integer_1d_array& p1,
     const ap::integer_1d_array& p2,
     const ap::real_1d_array& aoriginal,
     int n,
     bool& waserrors)
{
    int i;
    ap::real_1d_array a2;
    double t;
    ap::integer_1d_array f;

    a2.setbounds(0, n-1);
    f.setbounds(0, n-1);
    
    //
    // is set ordered?
    //
    for(i = 0; i <= n-2; i++)
    {
        waserrors = waserrors||ap::fp_greater(asorted(i),asorted(i+1));
    }
    
    //
    // P1 correctness
    //
    for(i = 0; i <= n-1; i++)
    {
        waserrors = waserrors||ap::fp_neq(asorted(i),aoriginal(p1(i)));
    }
    for(i = 0; i <= n-1; i++)
    {
        f(i) = 0;
    }
    for(i = 0; i <= n-1; i++)
    {
        f(p1(i)) = f(p1(i))+1;
    }
    for(i = 0; i <= n-1; i++)
    {
        waserrors = waserrors||f(i)!=1;
    }
    
    //
    // P2 correctness
    //
    for(i = 0; i <= n-1; i++)
    {
        a2(i) = aoriginal(i);
    }
    for(i = 0; i <= n-1; i++)
    {
        if( p2(i)!=i )
        {
            t = a2(i);
            a2(i) = a2(p2(i));
            a2(p2(i)) = t;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        waserrors = waserrors||ap::fp_neq(asorted(i),a2(i));
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testbdssunit_test_silent()
{
    bool result;

    result = testbdss(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testbdssunit_test()
{
    bool result;

    result = testbdss(false);
    return result;
}





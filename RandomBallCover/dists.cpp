/* This file is part of the Random Ball Cover (RBC) library.
 * (C) Copyright 2011, Lawrence Cayton [lcayton@tuebingen.mpg.de]
 */

#ifndef DISTS_C
#define DISTS_C

#include "dists.h"
#include "defs.h"

//computes the distance between the kth row of x and the lth row of y
real distVec(matrix x, matrix y, unint k, unint l){
  unint i, j;
  real sum=0;
  
  for(i=0; i<x.pc; i+=VEC_LEN){
    for(j=0; j<VEC_LEN; j++)
     sum += DIST( x.mat[IDX(k,i+j,x.ld)], y.mat[IDX(l,i+j,x.ld)] );
  }
  return DIST_EXP(sum);

}


//same as above, but will terminate early if the distance exceeds
//lb and return MAX_REAL
real distVecLB(matrix x, matrix y, unint k, unint l, real lb){
  unint i, j;
  real sum=0;
  real lb2 = lb==MAX_REAL? lb : DIST_ROOT(lb);

  for(i=0; i<x.pc; i+=VEC_LEN){
    for(j=0; j<VEC_LEN; j++)
      sum += DIST( x.mat[IDX(k,i+j,x.ld)], y.mat[IDX(l,i+j,x.ld)] );

    if( sum > lb2 ) 
      return MAX_REAL;

  }
  return DIST_EXP(sum);

}


#endif

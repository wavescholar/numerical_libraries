/* This file is part of the Random Ball Cover (RBC) library.
 * (C) Copyright 2011, Lawrence Cayton [lcayton@tuebingen.mpg.de]
 */

#ifndef BRUTE_C
#define BRUTE_C

#include "dists.h"
#include "utils.h"
#include "defs.h"
#include "brute.h"
#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include<gsl/gsl_sort.h>


// A basic parallel implementation of brute force 1-NN
void brutePar(matrix X, matrix Q, unint *NNs, real *dToNNs){
  real temp[CL];
  int i, j, k, t;
  
#pragma omp parallel for private(t,k,j,temp) 
  for( i=0; i<Q.pr/CL; i++ ){
    for(j=0; j<CL; j++)
      temp[j]=0; //hides compiler warning.
    
    t = i*CL;
    for(j=0; j<CL && t+j<Q.r;j++){
      dToNNs[t+j] = MAX_REAL;
      NNs[t+j] = 0;
    }
    
    for(j=0; j<X.r; j++ ){
      for(k=0; k<CL; k++){
	if(t+k<Q.r){
	  temp[k] = distVecLB( Q, X, t+k, j, dToNNs[t+k] );
	  //temp[k] = distVec( Q, X, t+k, j );
	}
      }
      for(k=0; k<CL; k++){
	if( t+k<Q.r && temp[k] < dToNNs[t+k]){
	  NNs[t+k] = j;
	  dToNNs[t+k] = temp[k];
	}
      }
    }
  }
}


// A basic implementation of brute force k-NN search.  This does not
// use a heap.  Instead it computes all distances and then sorts.
void bruteK(matrix x, matrix q, unint **NNs, real **dToNNs, unint k){
  int i, j, l;
  int nt = omp_get_max_threads();

  float ***d;
  size_t ***t;
  d = (float***)calloc(nt, sizeof(*d));
  t = (size_t***)calloc(nt, sizeof(*t));
  for(i=0; i<nt; i++){
    d[i] = (float**)calloc(CL, sizeof(**d));
    t[i] = (size_t**)calloc(CL, sizeof(**t));
    for(j=0; j<CL; j++){
      d[i][j] = (float*)calloc(x.pr, sizeof(***d));
      t[i][j] = (size_t*)calloc(x.pr, sizeof(***t));
    }
  }
      
#pragma omp parallel for private(j,l) shared(d,t,k)
  for( i=0; i<q.pr/CL; i++){
    int row = i*CL;
    int tn = omp_get_thread_num(); //thread num

    for( j=0; j<x.r; j++){
      for( l=0; l<CL; l++){
	d[tn][l][j] =  distVec( q, x, row+l, j);
      }
    }
    for(l=0; l<CL; l++)
      gsl_sort_float_smallest_index(t[tn][l], k, d[tn][l], 1, x.r);
    
    for(l=0; l<CL; l++){
      if(row+l<q.r){
	for(j=0; j<k; j++){
	  NNs[row+l][j] = (unint)t[tn][l][j];
	  dToNNs[row+l][j] = d[tn][l][t[tn][l][j]];
	}
      }
    }
  }

  for(i=0; i<nt; i++){
    for(j=0; j<CL; j++){
      free(d[i][j]);  free(t[i][j]);
    }
    free(d[i]);  free(t[i]);
  }  
  free(t); free(d);
}


// An implementation of brute k-NN search that uses a heap.
void bruteKHeap(matrix X, matrix Q, unint **NNs, real **dToNNs, unint K){
  real temp[CL];
  int i, j, k,t;

  int nt = omp_get_max_threads();
  heap **hp;
  hp = (heap**)calloc(nt, sizeof(*hp));
  for(i=0; i<nt; i++){
    hp[i] = (heap*)calloc(CL, sizeof(**hp));
    for(j=0; j<CL; j++)
      createHeap(&hp[i][j],K);
  }
  
#pragma omp parallel for private(t,k,j,temp) 
  for( i=0; i<Q.pr/CL; i++ ){
    t = i*CL;
    unint tn = omp_get_thread_num();
    heapEl newEl;

    for(j=0; j<X.r; j++ ){
      for(k=0; k<CL; k++){
	//temp[k] = distVec( Q, X, t+k, j );
	temp[k] = distVecLB( Q, X, t+k, j, hp[tn][k].h[0].val );
      }
      for(k=0; k<CL; k++){
	if( temp[k] <= hp[tn][k].h[0].val ){
	  newEl.id=j;
	  newEl.val=temp[k];
	  replaceMax( &hp[tn][k], newEl);
	}
      }
    }
    for(j=0; j<CL; j++)
      if(t+j<Q.r)
	heapSort(&hp[tn][j], NNs[t+j], dToNNs[t+j]);
      
    for(j=0; j<CL; j++)
      reInitHeap(&hp[tn][j]);
  }
  
  for(i=0; i<nt; i++){
    for(j=0; j<CL; j++)
      destroyHeap(&hp[i][j]);
    free(hp[i]);
  }
  free(hp);
}


// Performs a range count using brute force.
void rangeCount(matrix X, matrix Q, real *ranges, unint *counts){
  real temp;
  unint i, j, k;

#pragma omp parallel for private(j,k,temp) shared(counts,ranges)
  for( i=0; i<Q.pr/CL; i++ ){
    unint row = i*CL;
    for( j=0; j<CL; j++)
      counts[row+j] = 0;
    for(k=0; k<X.r; k++ ){
      for( j=0; j<CL; j++){
	temp = distVec( Q, X, row+j, k );
	counts[row+j] += ( temp < ranges[row+j] );
      }
    }
  }
}


// Performs a brute force NN search, but only between queries and points 
// belonging to each query's nearest representative.  This method is used by
// the one-shot algorithm.  
void bruteMap(matrix X, matrix Q, rep *ri, unint* qMap, unint *NNs, real *dToNNs){
  unint i, j, k;
  
  //Sort the queries, so that queries matched to a particular representative
  //will be processed together, improving cache performance.
  size_t *qSort = (size_t*)calloc(Q.pr, sizeof(*qSort));
  gsl_sort_uint_index(qSort,qMap,1,Q.r);

#pragma omp parallel for private(j,k) schedule(static)
  for( i=0; i<Q.pr/CL; i++ ){
    unint row = i*CL;
    for(j=0; j<CL; j++){
      dToNNs[qSort[row+j]] = MAX_REAL;
    }
    
    real temp;
    rep rt[CL];
    unint maxLen = 0;
    for(j=0; j<CL; j++){
      rt[j] = ri[qMap[qSort[row+j]]];
      maxLen = MAX(maxLen, rt[j].len);
    }  
    
    for(k=0; k<maxLen; k++ ){
      for(j=0; j<CL; j++ ){
	if( k<rt[j].len ){
	  temp = distVec( Q, X, qSort[row+j], rt[j].lr[k] ); //change to LB
	  if( temp < dToNNs[qSort[row+j]]){
	    NNs[qSort[row+j]] = rt[j].lr[k];
	    dToNNs[qSort[row+j]] = temp;
	  }
	}
      }
    }
  }
  free(qSort);
}


// Performs a brute force K-NN search, but only between queries and points 
// belonging to each query's nearest representative.  This method is used by
// the one-shot algorithm.  
void bruteMapK(matrix X, matrix Q, rep *ri, unint* qMap, unint **NNs, real **dToNNs, unint K){
  unint i, j, k;
  
  //Sort the queries, so that queries matched to a particular representative
  //will be processed together, improving cache performance.
  size_t *qSort = (size_t*)calloc(Q.pr, sizeof(*qSort));
  gsl_sort_uint_index(qSort,qMap,1,Q.r);

  unint nt = omp_get_max_threads();
  
  heap **hp;
  hp = (heap**)calloc(nt, sizeof(*hp));
  for(i=0; i<nt; i++){
    hp[i] = (heap*)calloc(CL, sizeof(**hp));
    for(j=0; j<CL; j++)
      createHeap(&hp[i][j],K);
  }   
  

#pragma omp parallel for private(j,k) schedule(static)
  for( i=0; i<Q.pr/CL; i++ ){
    unint row = i*CL;
    unint tn = omp_get_thread_num();
    heapEl newEl;
    
    real temp[CL];
    rep rt[CL];
    unint maxLen = 0;
    for(j=0; j<CL; j++){
      rt[j] = ri[qMap[qSort[row+j]]];
      maxLen = MAX(maxLen, rt[j].len);
      temp[j]=0.0; //gets rid of compiler warning
    }  
    
    for(j=0; j<maxLen; j++ ){
      for(k=0; k<CL; k++ ){
	if( j<rt[k].len )
	  temp[k] = distVecLB( Q, X, qSort[row+k], rt[k].lr[j], hp[tn][k].h[0].val  );
      }
      
      for(k=0; k<CL; k++ ){
	if( j<rt[k].len ){
	  if( temp[k] < hp[tn][k].h[0].val ){
	    newEl.id = rt[k].lr[j];
	    newEl.val = temp[k];
	    replaceMax( &hp[tn][k], newEl );
	  }
	}
      }
    }
    
    for(j=0; j<CL; j++)
      heapSort( &hp[tn][j], NNs[qSort[row+j]], dToNNs[qSort[row+j]] );
    
    for(j=0; j<CL; j++)
      reInitHeap(&hp[tn][j]);
  }

  for(i=0; i<nt; i++){
    for(j=0; j<CL; j++)
      destroyHeap(&hp[i][j]);
    free(hp[i]);
  }
  free(hp);
  free(qSort);
}


// Performs a brute force search between X and Q, but only compares distances
// between (x,q) pairs specified by toSearch.  This is used by the searchExact
// RBC method.
void bruteList(matrix X, matrix Q, rep *ri, intList *toSearch, unint numReps, unint *NNs, real *dToNNs){
  real temp;
  unint i, j, k, l;
  
  unint nt = omp_get_max_threads();
  unint m = Q.r;

  real **d = (real**)calloc(nt, sizeof(*d));
  unint **nn = (unint**)calloc(nt, sizeof(*nn));
  for(i=0; i<nt; i++){
    d[i] = (real*)calloc(m, sizeof(**d));
    nn[i] = (unint*)calloc(m, sizeof(**nn));
  }
  
  for(i=0; i<nt; i++){
    for(j=0; j<m; j++)
      d[i][j] = MAX_REAL;
  }


#pragma omp parallel for private(j,k,l,temp) //schedule(dynamic)
  for( i=0; i<numReps; i++ ){
    unint tn = omp_get_thread_num();

    for( j=0; j< toSearch[i].len/CL; j++){  //toSearch is assumed to be padded
      unint row = j*CL;
      unint qInd[CL];
      for(k=0; k<CL; k++)
	qInd[k] = toSearch[i].x[row+k];
      rep rt = ri[i];
      unint curMinInd[CL];
      real curMinDist[CL];
      for(k=0; k<CL; k++){
	curMinDist[k] = MAX_REAL;
	curMinInd[k] = 0; //hides compiler warning
      }
      for(k=0; k<rt.len; k++){
	for(l=0; l<CL; l++ ){
	  if(qInd[l]!=DUMMY_IDX){
	    temp = distVecLB( Q, X, qInd[l], rt.lr[k], curMinDist[l] );
	    if( temp <= curMinDist[l] ){
	      curMinInd[l] = rt.lr[k];
	      curMinDist[l] = temp;
	    }
	  }
	}
      }
      for(k=0; k<CL; k++){
	if(qInd[k]!=DUMMY_IDX && curMinDist[k] < d[tn][qInd[k]]){
	  nn[tn][qInd[k]] = curMinInd[k];
	  d[tn][qInd[k]] = curMinDist[k];
	}
      }
    }
  }

  // Could make this a parallel reduce....
  for( i=0; i<nt; i++){
    for( j=0; j<m; j++){
      if( d[i][j] <= dToNNs[j] ){
	dToNNs[j] = d[i][j];
	NNs[j] = nn[i][j];
      }
    }
  }
  for(i=0; i<nt; i++){
    free(d[i]); free(nn[i]);
  }
  free(d); free(nn);
}


// This method is the same as the above bruteList method, but for k-nn search.
// It uses a heap.
void bruteListK(matrix X, matrix Q, rep *ri, intList *toSearch, unint numReps, unint **NNs, real **dToNNs, unint K){
  real temp;
  unint i, j, k, l;
  
  unint nt = omp_get_max_threads();
  unint m = Q.r;

  heap **hp;
  hp = (heap**)calloc(nt, sizeof(*hp));
  for(i=0; i<nt; i++){
    hp[i] = (heap*)calloc(m, sizeof(**hp));
    for(j=0; j<m; j++)
      createHeap(&hp[i][j],K);
  }      
  
#pragma omp parallel for private(j,k,l,temp)
  for( i=0; i<numReps; i++ ){
    unint tn = omp_get_thread_num();
    heapEl newEl;
   
    for( j=0; j< toSearch[i].len/CL; j++){  //toSearch is assumed to be padded
      unint row = j*CL;
      unint qInd[CL];
      for(k=0; k<CL; k++)
	qInd[k] = toSearch[i].x[row+k];
      rep rt = ri[i];
      for(k=0; k<rt.len; k++){
	for(l=0; l<CL; l++ ){
	  if(qInd[l]!=DUMMY_IDX){
	    temp = distVecLB( Q, X, qInd[l], rt.lr[k], hp[tn][qInd[l]].h[0].val );
	    if( temp < hp[tn][qInd[l]].h[0].val ){
	      newEl.id = rt.lr[k];
	      newEl.val = temp;
	      replaceMax( &hp[tn][qInd[l]], newEl );
	    }
	  }
	}
      }
    }
  }
  
  // Now merge the NNs found by each thread.  Currently this performs the
  // merge within one thread, but should eventually be done in a proper
  // parallel-reduce fashion.
  unint **tInds;
  real **tVals;
  tInds = (unint**)calloc(nt, sizeof(*tInds));
  tVals = (real**)calloc(nt, sizeof(*tVals));
  for( i=0; i<nt; i++ ){
    tInds[i] = (unint*)calloc(K, sizeof(**tInds) );
    tVals[i] = (real*)calloc(K, sizeof(**tVals) );
  }
  
  size_t *indVec = (size_t*)calloc(nt*K, sizeof(*indVec));
  size_t *tempInds = (size_t*)calloc(nt*K, sizeof(*tempInds));
  real *valVec = (real*)calloc(nt*K, sizeof(*valVec));

  for( i=0; i<m; i++){
 
    for( j=0; j<nt; j++){
      heapSort( &hp[j][i], tInds[j], tVals[j] );
      for( k=0; k<K; k++){
	indVec[j*K + k] = tInds[j][k];
	valVec[j*K + k] = tVals[j][k];
      }
    }
    
    gsl_sort_float_index(tempInds, valVec, 1, nt*K);
    for( j=0; j<K; j++ ){
      dToNNs[i][j] = valVec[tempInds[j]];
      NNs[i][j] = indVec[tempInds[j]];
    }
  }

  free(tempInds);
  free(indVec);
  free(valVec);
  for( i=0; i<nt; i++ ){
    free(tInds[i]);
    free(tVals[i]);
  }
  free(tInds); 
  free(tVals);
  for(i=0; i<nt; i++){
    for(j=0; j<m; j++)
      destroyHeap(&hp[i][j]);
    free(hp[i]);
  }
  free(hp);
}
#endif

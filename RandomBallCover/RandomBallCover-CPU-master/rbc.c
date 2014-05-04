/* This file is part of the Random Ball Cover (RBC) library.
 * (C) Copyright 2011, Lawrence Cayton [lcayton@tuebingen.mpg.de]
 */

/* Core methods for building and searching the RBC. */

#ifndef RBC_C
#define RBC_C

#include "rbc.h"
#include "defs.h"
#include "utils.h"
#include "brute.h"
#include "dists.h"
#include<omp.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_sort.h>
#include<sys/time.h>
#include<stdio.h>
#include<math.h>
#include<stdint.h>


/* ************ EXACT SEARCH METHOD ************ */

//Builds the RBC for exact (1- or K-) NN search.
//Note: allocates memory for r and ri that must be freed
//externally.  Use freeRBC.
void buildExact(matrix x, matrix *r, rep *ri, unint numReps){
  unint n = x.r;
  unint i, j;
  unint longestLength = 0;

  if( numReps > n ){
    fprintf( stderr, "number of representatives must be less than the DB size\n");
    exit(1);
  }

  initMat(r, numReps, x.c);
  r->mat = (real*)calloc( sizeOfMat(*r), sizeof(*r->mat) );
  
  //pick r random reps
  pickReps(x,r);

  //Compute the rep for each x
  unint *repID = (unint*)calloc(x.pr, sizeof(*repID));
  real *dToReps = (real*)calloc(x.pr, sizeof(*dToReps));
  
  brutePar(*r,x,repID,dToReps);

  //gather the rep info & store it in struct
  for(i=0; i<numReps; i++){
    ri[i].len = 0;
    ri[i].radius = 0;
}    
  
  for(i=0; i<n; i++){
    ri[repID[i]].radius = MAX( dToReps[i], ri[repID[i]].radius );
    ri[repID[i]].len++;
  }
 
  unint **tempI = (unint**)calloc(numReps, sizeof(*tempI));
  real **tempD = (real**)calloc(numReps, sizeof(*tempD));
  for(i=0; i<numReps; i++){
    tempI[i] = (unint*)calloc(ri[i].len, sizeof(**tempI));
    tempD[i] = (real*)calloc(ri[i].len, sizeof(**tempD));
    ri[i].lr = (unint*)calloc(ri[i].len, sizeof(*ri[i].lr));
    ri[i].dists = (real*)calloc(ri[i].len, sizeof(*ri[i].dists));
    longestLength = MAX( longestLength, ri[i].len );
  }
  
  
  unint *tempCount = (unint*)calloc(numReps, sizeof(*tempCount));
  for(i=0; i<n; i++){
    tempI[repID[i]][tempCount[repID[i]]] = i;
    tempD[repID[i]][tempCount[repID[i]]++] = dToReps[i];
  }

  //this stores the owned points in order of distance to the
  //representative.  This ordering is not currently necc. for the search 
  //algorithm, but might be used in the future.
  size_t *p = (size_t*)calloc(longestLength, sizeof(*p));
  for(i=0; i<numReps; i++){
    gsl_sort_float_index( p, tempD[i], 1, ri[i].len );
    for(j=0; j<ri[i].len; j++){
      ri[i].dists[j] = tempD[i][p[j]];
      ri[i].lr[j] = tempI[i][p[j]];
    }
  }
  free(p);
  for(i=0;i<numReps; i++){
    free(tempI[i]);
    free(tempD[i]);
  }
  free(tempI);
  free(tempD);
  free(tempCount);
  free(dToReps);
  free(repID);
}



//Exact 1-NN search with the RBC.
void searchExact(matrix q, matrix x, matrix r, rep *ri, unint *NNs, real *dToNNs){
  unint i, j, k;
  unint *repID = (unint*)calloc(q.pr, sizeof(*repID));
  real *dToReps = (real*)calloc(q.pr, sizeof(*dToReps));
  intList *toSearch = (intList*)calloc(r.pr, sizeof(*toSearch));
  for(i=0;i<r.pr;i++)
    createList(&toSearch[i]);
  int nt = omp_get_max_threads();
  
  float ***d;  //d is indexed by: thread, cache line #, rep #
  d = (float***)calloc(nt, sizeof(*d));
  for(i=0; i<nt; i++){
    d[i] = (float**)calloc(CL, sizeof(**d));
    for(j=0; j<CL; j++){
      d[i][j] = (float*)calloc(r.pr, sizeof(***d));
    }
  }
  
#pragma omp parallel for private(j,k) //schedule(dynamic)
  for(i=0; i<q.pr/CL; i++){
    unint row = i*CL;
    unint tn = omp_get_thread_num();
    
    unint minID[CL];
    real minDist[CL];
    for(j=0;j<CL;j++)
      minDist[j] = MAX_REAL;
    
    for( j=0; j<r.r; j++ ){
      for(k=0; k<CL; k++){
	d[tn][k][j] = distVec(q, r, row+k, j);
	if(d[tn][k][j] < minDist[k]){
	  minDist[k] = d[tn][k][j]; //gamma
	  minID[k] = j;
	}
      }
    }

    for( j=0; j<CL; j++ )
      dToReps[row+j] = minDist[j];

    for(j=0; j<r.r; j++ ){
      for(k=0; k<CL; k++ ){
	real temp = d[tn][k][j];
	if( row + k<q.r && minDist[k] >= temp - ri[j].radius && 3.0*minDist[k] >= temp ){
#pragma omp critical
	  {
	    addToList(&toSearch[j], row+k);
	  }
	}
      }
    }
  }
  for(i=0; i<r.r; i++){
    while(toSearch[i].len % CL != 0)
      addToList(&toSearch[i],DUMMY_IDX);
  }

  bruteList(x,q,ri,toSearch,r.r,NNs,dToReps);
  
  for(i=0; i<q.r; i++)
    dToNNs[i] = dToReps[i];
  
  for(i=0;i<r.pr;i++)
    destroyList(&toSearch[i]);
  free(toSearch);
  free(repID);
  free(dToReps);
  for(i=0; i<nt; i++){
    for(j=0; j<CL; j++)
      free(d[i][j]); 
    free(d[i]);
  }
  free(d);
}


//Exact k-NN search with the RBC
void searchExactK(matrix q, matrix x, matrix r, rep *ri, unint **NNs, real **dNNs, unint K){
  unint i, j, k;
  unint *repID = (unint*)calloc(q.pr, sizeof(*repID));
  real **dToReps = (real**)calloc(q.pr, sizeof(*dToReps));
  for(i=0; i<q.pr; i++)
    dToReps[i] = (real*)calloc(K, sizeof(**dToReps));
  intList *toSearch = (intList*)calloc(r.pr, sizeof(*toSearch));
  for(i=0;i<r.pr;i++)
    createList(&toSearch[i]);
  int nt = omp_get_max_threads();

  float ***d;  //d is indexed by: thread, cache line #, rep #
  d = (float***)calloc(nt, sizeof(*d));
  for(i=0; i<nt; i++){
    d[i] = (float**)calloc(CL, sizeof(**d));
    for(j=0; j<CL; j++){
      d[i][j] = (float*)calloc(r.pr, sizeof(***d));
    }
  }
  
  heap **hp;
  hp = (heap**)calloc(nt, sizeof(*hp));
  for(i=0; i<nt; i++){
    hp[i] = (heap*)calloc(CL, sizeof(**hp));
    for(j=0; j<CL; j++)
      createHeap(&hp[i][j],K);
  }
  
#pragma omp parallel for private(j,k)
  for(i=0; i<q.pr/CL; i++){
    unint row = i*CL;
    unint tn = omp_get_thread_num();
    heapEl newEl; 

    for( j=0; j<r.r; j++ ){
      for(k=0; k<CL; k++){
	d[tn][k][j] = distVec(q, r, row+k, j);
	if( d[tn][k][j] < hp[tn][k].h[0].val ){
	  newEl.id = j;
	  newEl.val = d[tn][k][j];
	  replaceMax( &hp[tn][k], newEl );
	}
      }
    }
    for(j=0; j<r.r; j++ ){
      for(k=0; k<CL; k++ ){
	real minDist = hp[tn][k].h[0].val;
	real temp = d[tn][k][j];
	if( row + k<q.r && minDist >= temp - ri[j].radius && 3.0*minDist >= temp ){
#pragma omp critical
	  {
	    addToList(&toSearch[j], row+k);
	  }
	}
      }
    }
    for(j=0; j<CL; j++)
      reInitHeap(&hp[tn][j]);
  }

  for(i=0; i<r.r; i++){
    while(toSearch[i].len % CL != 0)
      addToList(&toSearch[i],DUMMY_IDX);
  }

  bruteListK(x,q,ri,toSearch,r.r,NNs,dNNs,K);

  
  //clean-up
  for(i=0; i<nt; i++){
    for(j=0; j<CL; j++)
      destroyHeap(&hp[i][j]);
    free(hp[i]);
  }
  free(hp);
  for(i=0;i<r.pr;i++)
    destroyList(&toSearch[i]);
  free(toSearch);
  free(repID);
  for(i=0;i<q.pr; i++)
    free(dToReps[i]);
  free(dToReps);
  for(i=0; i<nt; i++){
    for(j=0; j<CL; j++)
      free(d[i][j]); 
    free(d[i]);
  }
  free(d);
}


// Exact 1-NN search with the RBC.  This version works better on computers
// with a high core count (say > 4)
void searchExactManyCores(matrix q, matrix x, matrix r, rep *ri, unint *NNs, real *dToNNs){
  unint i, j, k;
  unint *repID = (unint*)calloc(q.pr, sizeof(*repID));
  real *dToReps = (real*)calloc(q.pr, sizeof(*dToReps));
  intList *toSearch = (intList*)calloc(r.pr, sizeof(*toSearch));
  
  for(i=0;i<r.pr;i++)
    createList(&toSearch[i]);

  brutePar(r,q,repID,dToReps);

#pragma omp parallel for private(j,k)
  for(i=0; i<r.pr/CL; i++){
    unint row = CL*i;
    real temp[CL];
    
    for(j=0; j<q.r; j++ ){
      for(k=0; k<CL; k++){
	temp[k] = distVec( q, r, j, row+k );
      }
      for(k=0; k<CL; k++){
	//dToRep[j] is current UB on dist to j's NN
	//temp - ri[i].radius is LB to dist belonging to rep i
	if( row+k<r.r && dToReps[j] >= temp[k] - ri[row+k].radius && 3.0*dToReps[j] >= temp[k] )
	  addToList(&toSearch[row+k], j); //need to search rep 
      }
    }
    for(j=0;j<CL;j++){
      if(row+j<r.r){
	while(toSearch[row+j].len % CL != 0)
	  addToList(&toSearch[row+j],DUMMY_IDX);	
      }
    }
  }

  //Most of the time is spent in this method
  bruteList(x,q,ri,toSearch,r.r,NNs,dToReps);
  
  for(i=0; i<q.r; i++)
    dToNNs[i] = dToReps[i];
  
  for(i=0;i<r.pr;i++)
    destroyList(&toSearch[i]);
  free(toSearch);
  free(repID);
  free(dToReps);
}


// Exact k-NN search with the RBC.  This version works better on computers
// with a high core count (say > 4)
void searchExactManyCoresK(matrix q, matrix x, matrix r, rep *ri, unint **NNs, real **dNNs, unint K){
  unint i, j, k;
  unint **repID = (unint**)calloc(q.pr, sizeof(*repID));
  for(i=0; i<q.pr; i++)
    repID[i] = (unint*)calloc(K, sizeof(**repID));
  real **dToReps = (real**)calloc(q.pr, sizeof(*dToReps));
  for(i=0; i<q.pr; i++)
    dToReps[i] = (real*)calloc(K, sizeof(**dToReps));
  intList *toSearch = (intList*)calloc(r.pr, sizeof(*toSearch));
  for(i=0;i<r.pr;i++)
    createList(&toSearch[i]);
  
  bruteKHeap(r,q,repID,dToReps,K);

#pragma omp parallel for private(j,k)
  for(i=0; i<r.pr/CL; i++){
    unint row = CL*i;
    real temp[CL];
    
    for(j=0; j<q.r; j++ ){
      for(k=0; k<CL; k++){
	temp[k] = distVec( q, r, j, row+k );
      }
      for(k=0; k<CL; k++){
	//dToRep[j] is current UB on dist to j's NN
	//temp - ri[i].radius is LB to dist belonging to rep i
	if( row+k<r.r && 3.0*dToReps[j][K-1] >= temp[k] && dToReps[j][K-1] >= temp[k] - ri[row+k].radius) 
	  addToList(&toSearch[row+k], j); //query j needs to search rep 
      }
    }
    for(j=0;j<CL;j++){
      if(row+j<r.r){
	while(toSearch[row+j].len % CL != 0)
	  addToList(&toSearch[row+j],DUMMY_IDX);	
      }
    }
  }

  bruteListK(x,q,ri,toSearch,r.r,NNs,dToReps,K);
  
  for(i=0; i<q.r; i++){
    for(j=0; j<K; j++)
      dNNs[i][j]=dToReps[i][j];
  }

  for(i=0;i<q.pr;i++)
    free(dToReps[i]);
  free(dToReps);
  for(i=0;i<r.pr;i++)
    destroyList(&toSearch[i]);
  free(toSearch);
  for(i=0;i<q.pr;i++)
    free(repID[i]);
  free(repID);
}



/* ************ ONE-SHOT METHOD ************ */

//Builds the RBC for the One-shot (inexact) method.
//Note: allocates memory for r and ri that must be freed
//externally.  Use freeRBC.
void buildOneShot(matrix x, matrix *r, rep *ri, unint numReps){
  unint s = numReps; //number of points per rep. Set equal to numReps
                     //as suggested by theory. 
  unint ps = CPAD(s);
  unint i, j;
  
  if( numReps > x.r ){
    fprintf( stderr, "number of representatives must be less than the DB size\n");
    exit(1);
  }

  initMat( r, numReps, x.c );
  r->mat = (real*)calloc( sizeOfMat(*r), sizeof(*r->mat));

  //pick r random reps
  pickReps(x,r); 

  //Compute the ownership lists
  unint **repID = (unint**)calloc(r->pr, sizeof(*repID));
  real **dToNNs = (real**)calloc(r->pr, sizeof(*dToNNs));
  for( i=0; i<r->pr; i++){
    repID[i] = (unint*)calloc(CPAD(ps), sizeof(**repID));
    dToNNs[i] = (real*)calloc(CPAD(ps), sizeof(**dToNNs));
  }

  //need to find the radius such that each rep contains s points
  bruteKHeap(x,*r,repID,dToNNs,s);
  
  for( i=0; i<r->pr; i++){
    ri[i].lr = (unint*)calloc(ps, sizeof(*ri[i].lr));
    ri[i].dists = (real*)calloc(ps, sizeof(*ri[i].dists));
    ri[i].len = s;
    for (j=0; j<s; j++){
      ri[i].lr[j] = repID[i][j];
    }
    //ri[i].radius = distVec( *r, x, i, ri[i].lr[s-1]);  //Not actually needed by one-shot alg
  }
  
  for( i=0; i<r->pr; i++){
    free(dToNNs[i]);
    free(repID[i]);
  }
  free(dToNNs);
  free(repID);
}


// Performs (approx) 1-NN search with the RBC One-shot algorithm.
void searchOneShot(matrix q, matrix x, matrix r, rep *ri, unint *NNs){
  unint *repID = (unint*)calloc(q.pr, sizeof(*repID));
  real *dToReps = (real*)calloc(q.pr, sizeof(*dToReps));
  
  // Determine which rep each query is closest to.
  brutePar(r,q,repID,dToReps);
    
  // Search that rep's ownership list.
  bruteMap(x,q,ri,repID,NNs,dToReps);
  
  free(repID);
  free(dToReps);
}


// Performs (approx) K-NN search with the RBC One-shot algorithm.
void searchOneShotK(matrix q, matrix x, matrix r, rep *ri, unint **NNs, real **dNNs, unint K){

  unint *repID = (unint*)calloc(q.pr, sizeof(*repID));
  real *dT = (real*)calloc(q.pr, sizeof(*dT));
  
  // Determine which rep each query is closest to.
  brutePar(r,q,repID,dT);
  
  // Search that rep's ownership list.
  bruteMapK(x,q,ri,repID,NNs,dNNs,K);
  
  free(repID);
  free(dT);
}

/* ************ END ONE-SHOT  ************ */


// Chooses representatives at random from x and stores them in r.
void pickReps(matrix x, matrix *r){
  unint n = x.r;
  unint i, j;

  unint *shuf = (unint*)calloc(n, sizeof(*shuf));
  for(i=0; i<n; i++)
    shuf[i]=i;


  //generate a random permutation of 1..n
  struct timeval tv;
  gettimeofday(&tv,NULL);
  gsl_rng * rng;
  const gsl_rng_type *rngT;
  
  gsl_rng_env_setup();
  rngT = gsl_rng_default;
  rng = gsl_rng_alloc(rngT);
  gsl_rng_set(rng,tv.tv_usec);
  
  gsl_ran_shuffle(rng, shuf, n, sizeof(*shuf));
  gsl_rng_free(rng);
 
  for(i=0; i<r->r; i++){
    for(j=0; j<r->c; j++){
      r->mat[IDX( i, j, r->ld )] = x.mat[IDX( shuf[i], j, x.ld )];
    }
  }
  free(shuf);
}


// Determines the total number of computations needed by RBC to 
// find the the NNS.  This function is useful mainly for 
// evaluating the effectiveness of the RBC.  
void searchStats(matrix q, matrix x, matrix r, rep *ri, double *avgDists){
  unint i, j;
  unint *repID = (unint*)calloc(q.pr, sizeof(*repID));
  real *dToReps = (real*)calloc(q.pr, sizeof(*dToReps));

  brutePar(r,q,repID,dToReps);

  //for each q, need to determine which reps to examine
  size_t numAdded=0;
  size_t totalComp=0;
#pragma omp parallel for private(j) reduction(+:numAdded,totalComp)
  for(i=0; i<q.r; i++){
    for(j=0; j<r.r; j++ ){
      real temp = distVec( q, r, i, j );
      //dToRep[i] is current UB on dist to i's NN
      //temp - ri[j].radius is LB to dist belonging to rep j
      if( dToReps[i] >= temp - ri[j].radius && temp <= 3.0*dToReps[i] ){
	numAdded++;
	totalComp+=ri[j].len;
      }
    }
  }
  
  *avgDists = ((double)totalComp)/((double)q.r);
  free(repID);
  free(dToReps);
}


//frees the memory associated with the RBC
void freeRBC(matrix r, rep *ri){
  unint i;

  free(r.mat);
  for( i=0; i<r.pr; i++ ){
    free( ri[i].lr );
    if( ri[i].dists )
      free( ri[i].dists );
  }
  free(ri);

}

#endif

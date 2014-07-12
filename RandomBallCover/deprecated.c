//This file collects versions of functions that are no longer used.

//Returns a length l subset of a random permutation of [0,...,n-1]
//using the knuth shuffle.
//Input variable x is assumed to be alloced and of size l.

void subRandPerm(unint,unint,unint*);
void subRandPerm(unint l, unint n, unint *x){
  unint i,ri, *y;
  y = (unint*)calloc(n,sizeof(*y));
    
  struct timeval t3;
  gettimeofday(&t3,NULL);
  srand(t3.tv_usec);
  
  for(i=0;i<n;i++)
    y[i]=i;
  
  for(i=0;i<MIN(l,n-1);i++){  //The n-1 bit is necessary because you can't swap the last 
                              //element with something larger.
    ri=randBetween(i+1,n);
    swap(&y[i],&y[ri]);
  }
  
  for(i=0;i<l;i++)
    x[i]=y[i];
  free(y);
}


void bruteList(matrix,matrix,rep*,intList*,unint,unint*,real*);
void bruteList(matrix X, matrix Q, rep *ri, intList *toSearch, unint numReps, unint *NNs, real *dToNNs){
  real temp;
  unint i, j, k, l;
  
  for(i=0; i<Q.r; i++){
    dToNNs[i] = MAX_REAL;
    NNs[i] = 0;
  }

#pragma omp parallel for private(j,k,l,temp)
  for( i=0; i<numReps; i++ ){
    for( j=0; j< toSearch[i].len/CL; j++){  //toSearch is assumed to be padded
      unint row = j*CL;
      unint qInd[CL];
      for(k=0; k<CL; k++)
	qInd[k] = toSearch[i].x[row+k];
      rep rt = ri[i];
      unint curMinInd[CL];
      real curMinDist[CL];
      for(k=0; k<CL; k++)
	curMinDist[k] = MAX_REAL;
      for(k=0; k<rt.len; k++){
	for(l=0; l<CL; l++ ){
	  if(qInd[l]!=DUMMY_IDX){
	    temp = distVec( Q, X, qInd[l], rt.lr[k] );
	    if( temp < curMinDist[l] ){
	      curMinInd[l] = rt.lr[k];
	      curMinDist[l] = temp;
	    }
	  }
	}
      }
#pragma omp critical
      {
	for(k=0; k<CL; k++){
	  if( qInd[k]!=DUMMY_IDX && curMinDist[k] < dToNNs[qInd[k]]){
	    NNs[qInd[k]] = curMinInd[k];
	    dToNNs[qInd[k]] = curMinDist[k];
	  }
	}
      }
    }
  }
}


void bruteKDists(matrix,matrix,size_t**,real**,unint);
void bruteKDists(matrix x, matrix q, size_t **NNs, real **D, unint k){
  int i, j;

  float **d;
  d = (float**)calloc(q.pr, sizeof(*d));
  size_t **t;
  t = (size_t**)calloc(q.pr, sizeof(*t));
  for( i=0; i<q.pr; i++){
    d[i] = (float*)calloc(x.pr, sizeof(**d));
    t[i] = (size_t*)calloc(x.pr, sizeof(**t));
  }

#pragma omp parallel for private(j)
  for( i=0; i<q.r; i++){
    for( j=0; j<x.r; j++)
      d[i][j] = distVec( q, x, i, j );
    gsl_sort_float_index(t[i], d[i], 1, x.r);
    for ( j=0; j<k; j++){
      NNs[i][j]=t[i][j];
      D[i][j]=d[i][t[i][j]];
    }
  }

  for( i=0; i<q.pr; i++){
    free(t[i]);
    free(d[i]);
  }
  free(t);
  free(d);
}

void rangeCount(matrix,matrix,real*,unint*);
void rangeCount(matrix X, matrix Q, real *ranges, unint *counts){
  real temp;
  unint i, j;

#pragma omp parallel for private(j,temp)
  for( i=0; i<Q.r; i++ ){
    counts[i] = 0;
    for(j=0; j<X.r; j++ ){
      temp = distVec( Q, X, i, j );
      counts[i] += ( temp < ranges[i] );
    }
  }
}


void bruteMap(matrix,matrix,rep*,unint*,unint*,real*);
//This is a very crude implementation without any reordering of q
void bruteMap(matrix X, matrix Q, rep *ri, unint* qMap, unint *NNs, real *dToNNs){
  real temp;
  unint i, j;
  
#pragma omp parallel for private(j,temp)
  for( i=0; i<Q.r; i++ ){
    dToNNs[i] = MAX_REAL;
    NNs[i] = 0;
    rep rt = ri[qMap[i]];
    for(j=0; j<rt.len; j++ ){
      temp = distVec( Q, X, i, rt.lr[j] );
      if( temp < dToNNs[i]){
	NNs[i] = rt.lr[j];
	dToNNs[i] = temp;
      }
    }
  }
}

void bruteMap2(matrix,matrix,rep*,unint*,unint*,real*);
//This is a still crude, but doesn't ignore the cache at least
void bruteMap2(matrix X, matrix Q, rep *ri, unint* qMap, unint *NNs, real *dToNNs){
  unint i, j, k;
  
#pragma omp parallel for private(j,k)
  for( i=0; i<Q.pr/CL; i++ ){
    unint row = i*CL;
    for(j=0; j<CL; j++){
      dToNNs[row+j] = MAX_REAL;
      NNs[row+j] = 0;
    }
    
    real temp;
    rep rt[CL];
    unint maxLen = 0;
    for(j=0; j<CL; j++){
      rt[j] = ri[qMap[row+j]];
      maxLen = MAX(maxLen, rt[j].len);
    }  
    
    for(k=0; k<maxLen; k++ ){
      for(j=0; j<CL; j++ ){
	if( k<rt[j].len ){
	  temp = distVec( Q, X, row+j, rt[j].lr[k] );
	  if( temp < dToNNs[row+j]){
	    NNs[row+j] = rt[j].lr[k];
	    dToNNs[row+j] = temp;
	  }
	}
      }
    }
  }
}

void bruteK(matrix,matrix,size_t**,unint);
//totally unoptimized
void bruteK(matrix x, matrix q, size_t **NNs, unint k){
  int i, j;

  float **d;
  d = (float**)calloc(q.pr, sizeof(*d));
  size_t **t;
  t = (size_t**)calloc(q.pr, sizeof(*t));
  for( i=0; i<q.pr; i++){
    d[i] = (float*)calloc(x.pr, sizeof(**d));
    t[i] = (size_t*)calloc(x.pr, sizeof(**t));
  }

#pragma omp parallel for private(j)
  for( i=0; i<q.r; i++){
    for( j=0; j<x.r; j++)
      d[i][j] = distVec( q, x, i, j );
    gsl_sort_float_index(t[i], d[i], 1, x.r);
    for ( j=0; j<k; j++)
      NNs[i][j]=t[i][j];
  }

  for( i=0; i<q.pr; i++){
    free(t[i]);
    free(d[i]);
  }
  free(t);
  free(d);
}
void brute(matrix,matrix,unint*,real*);
void brute(matrix X, matrix Q, unint *NNs, real *dToNNs){
  real temp;
  unint i, j;

  for( i=0; i<Q.r; i++ ){
    dToNNs[i] = MAX_REAL;
    NNs[i] = 0;
    for(j=0; j<X.r; j++ ){
      temp = distVec( Q, X, i, j );
      if( temp < dToNNs[i]){
	NNs[i] = j;
	dToNNs[i] = temp;
      }
    }
  }
}

void bruteCache(matrix,matrix,unint*,real*);
void bruteCache(matrix X, matrix Q, unint *NNs, real *dToNNs){
  real temp[CL];

  int i, j, k, t;
  for( i=0; i<Q.pr/CL; i++ ){
    t = i*CL;
    for(j=0;j<CL;j++){
      dToNNs[t+j] = MAX_REAL;
      NNs[t+j] = 0;
    }
    for(j=0; j<X.r; j++ ){
      for(k=0; k<CL; k++){
	temp[k] = distVec( Q, X, t+k, j );
      }
      for(k=0; k<CL; k++){
	if( temp[k] < dToNNs[t+k]){
	  NNs[t+k] = j;
	  dToNNs[t+k] = temp[k];
	}
      }
    }
  }
}


void brutePar(matrix,matrix,unint*,real*);
// Most basic implementation of parallel brute force search.  Not very fast.
void brutePar(matrix X, matrix Q, unint *NNs, real *dToNNs){
  real temp;
  int i, j;
 
  
#pragma omp parallel for private(j,temp)
  for( i=0; i<Q.r; i++ ){
    dToNNs[i] = MAX_REAL;
    NNs[i] = 0;
    for(j=0; j<X.r; j++ ){
      temp = distVec( Q, X, i, j );
      if( temp < dToNNs[i]){
	NNs[i] = j;
	dToNNs[i] = temp;
      }
    }
  }
}



void bruteList(matrix,matrix,rep*,intList*,unint*,real*);
void bruteList(matrix X, matrix Q, rep *ri, intList *toSearch, unint *NNs, real *dToNNs){
  real temp;
  unint i, j, k;
  
#pragma omp parallel for private(j,k,temp)
  for( i=0; i<Q.r; i++ ){
    dToNNs[i] = MAX_REAL;
    NNs[i] = 0;
    for( j=0; j<toSearch[i].len; j++ ){
      rep rt = ri[ toSearch[i].x[j] ];
      for(k=0; k<rt.len; k++ ){
	temp = distVec( Q, X, i, rt.lr[k] );
	if( temp < dToNNs[i]){
	  NNs[i] = rt.lr[k];
	  dToNNs[i] = temp;
	}
      }
    }
  }
}


void dimEst(matrix,unint);
void dimEst(matrix x, unint s){
  unint n = x.r;
  unint i,j ;
  int t = 3;
  long unsigned int curMax, oldMax; 

  matrix r;
  r.c=x.c; r.pc=x.pc; r.r=s; r.pr=CPAD(s); r.ld=r.pc;
  r.mat = (real*)calloc( r.pc*r.pr, sizeof(*r.mat) );
  
  size_t** NNs = (size_t**)calloc(r.r, sizeof(*NNs));
  real** dToNNs = (real**)calloc(r.r, sizeof(*dToNNs));
  for(i=0;i<r.pr;i++){
    NNs[i] = (size_t*)calloc(t, sizeof(**NNs));
    dToNNs[i] = (real*)calloc(t, sizeof(**dToNNs));
  }
  
  real* ranges = (real*)calloc(r.pr, sizeof(*ranges));
  unint *counts = (unint*)calloc(r.pr,sizeof(*counts));

  size_t *p = (size_t*)calloc(r.pr, sizeof(*p));
  //pick r random reps
  pickReps(x,&r);

  //printf("calling bruteKdists\n");
  bruteKDists(x,r,NNs,dToNNs,t); 
  //  printf("done\n");
  for(i=0;i<r.r;i++)
    ranges[i] = dToNNs[i][t-1];

  for(i=0; i<10; i++){
    rangeCount(x,r,ranges,counts);
 
    //gsl_sort_uint_index(p,counts,1,r.r);
    //printf("80 = %u \n",counts[p[5*r.r/10]]);
    for(j=0; j<r.r; j++)
      ranges[j]*=2.0;
    curMax = 0;
    unint avg = 0;
    for(j=0; j<r.r; j++){
      //      printf("%u ",counts[j]);
      curMax = MAX(counts[j],curMax);
      avg +=counts[j];
    }
    //    printf("\n");
    //    printf("avg = %6.4f \n",((double)avg)/((double)r.r));
    //    printf("%lu \n",curMax);
  }
  
  for(i=0;i<r.r;i++){
    free(NNs[i]);
    free(dToNNs[i]);
  }
  free(NNs);
  free(dToNNs);
  free(r.mat);
  free(ranges);
  free(counts);
  free(p);
}


void searchExact2(matrix,matrix,matrix,rep*,unint*);

void searchExact2(matrix q, matrix x, matrix r, rep *ri, unint *NNs){
  unint i, j;
  unint *repID = (unint*)calloc(q.pr, sizeof(*repID));
  real *dToReps = (real*)calloc(q.pr, sizeof(*dToReps));
  intList *toSearch = (intList*)calloc(r.pr, sizeof(*toSearch));
  for(i=0;i<r.pr;i++)
    createList(&toSearch[i]);
  
  struct timeval tvB,tvE;

  gettimeofday(&tvB,NULL);
  brutePar2(r,q,repID,dToReps);
  gettimeofday(&tvE,NULL);
  printf("[SE]brutePar2 time elapsed = %6.4f \n", timeDiff(tvB,tvE) );

  gettimeofday(&tvB,NULL);
#pragma omp parallel for private(j)
  for(i=0; i<r.r; i++){
    for(j=0; j<q.r; j++ ){
      real temp = distVec( q, r, j, i );
      //dToRep[j] is current UB on dist to j's NN
      //temp - ri[i].radius is LB to dist belonging to rep i
      if( dToReps[j] >= temp - ri[i].radius)
	addToList(&toSearch[i], j); //need to search rep i
    }
    while(toSearch[i].len % CL != 0)
      addToList(&toSearch[i],DUMMY_IDX);
  }
  gettimeofday(&tvE,NULL);
  printf("[SE]loop time elapsed = %6.4f \n", timeDiff(tvB,tvE) );

  gettimeofday(&tvB,NULL);
  bruteListRev(x,q,ri,toSearch,r.r,NNs,dToReps);
  gettimeofday(&tvE,NULL);
  printf("[SE]bruteListRev time elapsed = %6.4f \n", timeDiff(tvB,tvE) );

  for(i=0;i<r.pr;i++)
    destroyList(&toSearch[i]);
  free(toSearch);
  free(repID);
  free(dToReps);
}

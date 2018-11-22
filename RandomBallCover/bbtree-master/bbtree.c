#ifndef BBTREE_C
#define BBTREE_C

//This file is part of the bregman ball tree package.
//(c) copyright 2008, Lawrence Cayton

#include<limits.h>
#include<stdlib.h>
#include<stdio.h>
#include<sys/time.h>
#include<math.h>
#include "bbtree.h"
#include "utils.h"

treenode* buildTree(double **data, int n, int d, bregdiv div,int bucketSize){
  int i;
  /*seed the random number generater */
  struct timeval t3;
  gettimeofday(&t3,NULL);
  srand(t3.tv_usec);
  
  int* inds;
  inds= calloc(n,sizeof(int));
  
  if(inds==NULL)
    printf("error with allocation");
  for(i=0;i<n;i++)
    inds[i]=i;

  treenode* root = rBuild(data,inds,n,d,div,bucketSize);
  return root;
  free(inds);
}

treenode* rBuild(double **data, int* inds, int n, int d, bregdiv div, int bucketSize){
  int i;
  int *indsL, *indsR;
  double divTemp;

  indsL = calloc(n,sizeof(int));
  indsR = calloc(n,sizeof(int)); 
  
  treenode* node = setupNode(n,d);
  node->n=n;
  for(i=0;i<n;i++)
    node->inds[i]=inds[i];
  mean(node->mu,data,inds,n,d);
  
  node->R=0;
  for(i=0;i<n;i++){
    divTemp = div.div(data[inds[i]],node->mu,d);
    if(divTemp> node->R)
      node->R=divTemp;
  }
      
  div.gradf(node->mup,node->mu,d); 
  
  if (n<= bucketSize){
    node->isLeaf = 1;
    return node;
  }
  node->isLeaf=0;
  int lLength=0, rLength=0,its=0;
  
  while(its<10 && (lLength==0 || rLength==0)){
    findSplitKmeans(indsL,indsR,&lLength,&rLength,data,inds,n,d,div);
    its++;
  }
  
  //findSplit(indsL,indsR,&lLength,&rLength,data,inds,n,d,div);
  //findSplitMaxVar(indsL,indsR,&lLength,&rLength,data,inds,n,d,div);
  if(its < 10){
    node->l = rBuild(data,indsL,lLength,d,div,bucketSize);
    node->r = rBuild(data,indsR,rLength,d,div,bucketSize);
  }
  else{
    node->isLeaf=1;
    return node;
  }

  free(indsL);
  free(indsR);
  return node;
}

treenode* setupNode(int n, int d){
  treenode *node;
  node = (treenode*)calloc(1,sizeof(treenode));
  node->mu = calloc(d,sizeof(double));
  node->mup = calloc(d,sizeof(double));
  node->inds = calloc(n,sizeof(int));
  
  return node;
}



void findSplitKmeans(int* indsL, int* indsR, int* lLength, int* rLength,double**data, int*inds, int n, int d, bregdiv div){
  int lp,rp;
  int i, rCount,lCount,j;
  double cost;

  int MAXITS = 10;
  double *mul, *mur;
  
  mul = calloc(d,sizeof(double));
  mur = calloc(d,sizeof(double));
  
  double *Dl,*Dr,*diff;
  Dl = calloc(n,sizeof(double));
  Dr = calloc(n,sizeof(double));
  diff = calloc(n,sizeof(double));


  lp = genRand(n); 
  while (lp == (rp = genRand(n))); 

  for(i=0;i<d;i++){ 
    mul[i] = data[inds[lp]][i]; 
    mur[i] = data[inds[rp]][i]; 
  } 
  
  for(j=0;j<MAXITS;j++){
    cost=0;
    calcDvectorNonInd(Dl,data,inds,n,d,mul,div);
    calcDvectorNonInd(Dr,data,inds,n,d,mur,div);
    lCount = 0;
    rCount = 0;
    for(i=0;i<n;i++){
      if(Dl[i] < Dr[i]){
	indsL[lCount++] = inds[i];
	cost+=Dl[i];
      }
      else{
	indsR[rCount++] = inds[i];
	cost+=Dr[i];
      }
    }
    if(lCount==0 || rCount==0){ //start over;
      //j=0;
      lp = genRand(n);
      while (lp == (rp = genRand(n)));
      for(i=0;i<d;i++){
	mul[i] = data[inds[lp]][i];
	mur[i] = data[inds[rp]][i];
      }
    }
    else{
      mean(mul,data,indsL,lCount,d); //On the last iteration, these mus
      mean(mur,data,indsR,rCount,d); //aren't used.
    }
  }
  
  (*lLength) = lCount;
  (*rLength) = rCount;
    
  free(mul);
  free(mur);
  free(diff);
  free(Dl);
  free(Dr);

}

void findSplit(int* indsL, int* indsR, int* lLength, int* rLength,double**data, int*inds, int n, int d, bregdiv div){
  int p,lp,rp;
  int i;

  *lLength = floor(n/2.0);
  *rLength = ceil(n/2.0);

  p = genRand(n);
  lp = calcFarthest(data,inds,n,d,p,div);
  rp = calcFarthest(data,inds,n,d,lp,div);
    
  double *Dl,*Dr,*diff;
  Dl = calloc(n,sizeof(double));
  Dr = calloc(n,sizeof(double));
  diff = calloc(n,sizeof(double));
  
  calcDvector(Dl,data,inds,n,d,lp,div);
  calcDvector(Dr,data,inds,n,d,rp,div);
  for(i=0;i<n;i++)
    diff[i]=Dl[i]-Dr[i];
  
  sortWithInds(diff,inds,n);
  
  for(i=0;i<*lLength;i++){
    indsL[i]=inds[i];
  }
  

  for(i=0;i<*rLength;i++){
    indsR[i]=inds[i+*lLength];
  }

  free(diff);
  free(Dl);
  free(Dr);

}

void findSplitMaxVar(int* indsL, int* indsR, int* lLength, int* rLength,double**data, int*inds, int n, int d, bregdiv div){
  int rp;
  int i,j;

  *lLength = floor(n/2.0);
  *rLength = ceil(n/2.0);

  //  p = genRand(n);
  //  lp = calcFarthest(data,inds,n,d,p,div);
  //  rp = calcFarthest(data,inds,n,d,lp,div);
    
  double *Dl,*Dr,*diff,var;
  Dl = calloc(n,sizeof(double));
  Dr = calloc(n,sizeof(double));
  diff = calloc(n,sizeof(double));
  double maxSoFar=0;
  int bestLP=0;
  for(i=0;i<n;i++){
    rp = calcFarthest(data,inds,n,d,i,div);
    var =0;
    calcDvector(Dl,data,inds,n,d,i,div);
    calcDvector(Dr,data,inds,n,d,rp,div);
    for(j=0;j<n;j++)
      var+=fabs(Dl[j]-Dr[j]);
    if(var>maxSoFar){
      maxSoFar=var;
      bestLP=i;
    }
  }
  
  rp = calcFarthest(data,inds,n,d,bestLP,div);
  calcDvector(Dl,data,inds,n,d,bestLP,div);
  calcDvector(Dr,data,inds,n,d,rp,div);
  for(i=0;i<n;i++)
    diff[i]=Dl[i]-Dr[i];
  sortWithInds(diff,inds,n);
  
  for(i=0;i<*lLength;i++){
    indsL[i]=inds[i];
  }
  
  for(i=0;i<*rLength;i++){
    indsR[i]=inds[i+*lLength];
  }

  free(diff);
  free(Dl);
  free(Dr);

}



// This helper returns the index of the farthest pt to pt p
int calcFarthest(double **data, int*inds, int n, int d, int p, bregdiv div){
  int i;
  int curBest=0;
  double curDiv=0,tempDiv;
  
  for(i=0;i<n;i++){
    tempDiv = div.div(data[inds[i]],data[inds[p]],d);
    if (tempDiv > curDiv){
      curDiv=tempDiv;
      curBest=i;
    }
  }
  return curBest;
}

// This is a helper function that calculates the 
// divergence from a set of pts to point p, where p is an index
void calcDvector(double *D, double **data, int* inds,int n,int d,int p, bregdiv div){
  int i;
  for(i=0;i<n;i++)
    D[i] = div.div(data[inds[i]],data[inds[p]],d);
}

void calcDvectorNonInd(double *D, double **data, int* inds,int n,int d,double *p, bregdiv div){
  int i;
  for(i=0;i<n;i++)
    D[i] = div.div(data[inds[i]],p,d);
}

void deleteTree(treenode* root){
  
  if(!(root->isLeaf)){
    deleteTree(root->l);
    deleteTree(root->r);
  }
  
  free(root->mu);
  free(root->mup);
  free(root->inds);
  free(root);

}

#endif


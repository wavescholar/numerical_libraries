#ifndef SEARCH_C
#define SEARCH_C

//This file is part of the bregman ball tree package.
//(c) copyright 2008, Lawrence Cayton

#include "bregdiv.h"
#include "bbtree.h"
#include "search.h"
#include "utils.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

double *x, *xp;  //temp variables used throughout
long distCount;
double *dtoNNs;
int *NNs;
int leafBound=INT_MAX;


//This method performs nearest neighbor search for m queries.
void multisearch(treenode* root, double **Q, double**data, bregdiv div, int n, int d, int m, double eps, int k, int maxLeaves){
  int i,j;
  double* qp;
  x = calloc(d,sizeof(double));
  xp = calloc(d,sizeof(double));
  qp = calloc(d,sizeof(double));
  dtoNNs = calloc(k,sizeof(double));
  NNs = calloc(k,sizeof(int));

  //Note that dtoNNs and NNs should really be implemented using a 
  //efficient heap structure.  They are currently just kept as sorted
  //arrays, which is probably fine if k is small.

  distCount=0;
  for(i=0;i<m;i++){
    leafBound = maxLeaves;
    for(j=0;j<k;j++){
      dtoNNs[j]=HUGE_VAL;
      NNs[j]=-1;
    }
        
    div.gradf(qp,Q[i],d);
    rsearch(root,Q[i],qp,data,div,d,eps,k);
    
    /*printf("query %d nns are \n",i);
    for(j=0;j<k;j++)
      printf("%d ",NNs[j]);
      printf("\n"); */
  }

  free(x);
  free(xp);
  free(qp);
  free(NNs);
  free(dtoNNs);

}

void rsearch(treenode *node, double *q, double *qp, double**data, bregdiv div, int d, double eps, int k){
  double divTemp;
  int i;

  if(node->isLeaf){ //calc dists to all pts.
    distCount+=node->n;
    for(i=0;i<(node->n);i++){
      divTemp = div.div(data[node->inds[i]],q,d);
      if (NNs[0]==-1 || divTemp< dtoNNs[0]){
	//put the new point on the "heap"
	//	printf("inserting %d \n",node->inds[i]);
	insert(NNs,dtoNNs,node->inds[i],divTemp,k);  
	
      } 
    }
    return;
  }
  
  distCount+=2;
  node->l->divToMu = div.div(q,node->l->mu,d);
  node->r->divToMu = div.div(q,node->r->mu,d);
  if(node->l->divToMu < node->r->divToMu){ //search left
    rsearch(node->l,q,qp,data,div,d,eps,k);
    if(needToSearch(node->r,q,qp,div,d,dtoNNs[0],eps))
      rsearch(node->r,q,qp,data,div,d,eps,k);
  }
  else{ //search right
    rsearch(node->r,q,qp,data,div,d,eps,k);
    if(needToSearch(node->l,q,qp,div,d,dtoNNs[0],eps))
      rsearch(node->l,q,qp,data,div,d,eps,k);
  }
}


int needToSearch(treenode *node,double* q,double* qp,bregdiv div, int d, double minDistEst,double eps){
  /*
    //In l2^2 case, can use the following (derived from triangle inequality)

    if( node->divToMu < node->R)
      return 1;
    else if(minDistEst < 0.5*(sqrt(2.0*node->divToMu) - sqrt(2.0*node->R))*(sqrt(2.0*node->divToMu) - sqrt(2.0*node->R)))
      return 0;
    else
      return 1; 
  */

    if( ((node->divToMu) < (node->R)) || ((node->divToMu)*(1.0+eps) < (minDistEst)))
      return 1;
    return recBinSearch(0.0,1.0,node,q,qp,div,d,minDistEst,eps);
}

int recBinSearch(double l, double r, treenode *node, double *q, double* qp, bregdiv div, int d, double minDistEst,double eps){
  double theta = 0.5*(double)(l+r);
  double divToMu, divToQ, lbnd;
  int i;

  distCount+=2;

  for(i=0;i<d;i++){
    xp[i]=theta*qp[i]+(1-theta)*(node->mup[i]);
  }
  div.gradfstar(x,xp,d);
  
  //evaluate dual for new lower bnd
  divToMu = div.div(x,node->mu,d);
  divToQ = div.div(x,q,d);
  lbnd = divToQ + (1.0/theta -1.0)*(divToMu - node->R);
  
  if((1.0+eps)*lbnd >= minDistEst)
    return 0; 
    
  if( (fabs(divToMu - node->R)) < node->R*CLOSEENOUGH) 
      return 1;
  
  if(divToMu > node->R)
    return recBinSearch(l,theta,node,q,qp,div,d,minDistEst,eps);
  else{ //divToQ is an upper bnd, so we might be able to stop
    if(divToQ*(1.0+eps) < minDistEst){
      return 1;
    }
    return recBinSearch(theta,r,node,q,qp,div,d,minDistEst,eps);
  }
}


#endif

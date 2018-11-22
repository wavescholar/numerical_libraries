#ifndef UTILS_C
#define UTILS_C

//This file is part of the bregman ball tree package.
//(c) copyright 2008, Lawrence Cayton

#include<stdarg.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include "utils.h"
#include "bbtree.h"

void mean(double* mu,double** data,int* inds, int n, int d){
  int i,j;
  
  for(i=0;i<d;i++){
    for(j=0;j<n;j++){
      mu[i]+=data[inds[j]][i];
    }
    mu[i]/=((double)n);
  }
}


/* generates a random int in [0,n-1] */
/* note: assumes that rand has already been seeded */
int genRand(int n){
  int val;
  while(n<= (val= rand()/(int)(((unsigned)RAND_MAX + 1) / n)));
  return val;
}


/* ****************************** */
/* the following functions adapt qsort to return the sorted indices */
int compInds(const void *a, const void *b){
  double diff = (double)( (*((sortitem*)a)).val - (*((sortitem*)b)).val );
  if(diff<0)
    return -1;
  else if(diff==0)
    return 0;
  else
    return 1;
}

void sortWithInds(double* a, int* inds, int length){
  int i;
  sortitem *s;
  s = calloc(length,sizeof(sortitem));

  for(i=0;i<length;i++){
    s[i].val = a[i];
    s[i].ind = inds[i];
  }

  qsort(s,length,sizeof(sortitem),compInds);
  
  for(i=0;i<length;i++){
    inds[i]=s[i].ind;
    a[i]=s[i].val;
  }
    
  free(s);
}
/* ****************************** */

void readData(double **x, int n, int d, char* file){
  int i,j;
  float t;
  FILE *fp = fopen(file,"r");
  
  if(fp==NULL){
    fprintf(stderr,"error opening file.. exiting\n");
    exit(1);
  }

  for(i=0;i<n;i++){
    for(j=0;j<d;j++){
      if(fscanf(fp,"%f ", &t)==EOF){
	fprintf(stderr,"error reading file.. exiting \n");
	exit(1);
      }
      x[i][j]=(double)t;
    }
  }
  fclose(fp);
}


void writeDoubs(int num, char* file, double x,...){
  double z;
  int i;

  FILE*fp = fopen(file,"a");
  va_list s;
  va_start(s,x);
  
  for(z=x, i=0; i<num; z=va_arg(s,double),i++)
    fprintf(fp,"%6.5f ",z);
  fprintf(fp,"\n");
  fclose(fp);
}


void writeTree(treenode* root,int d, char* file){
  FILE *fp = fopen(file,"w");
  if(fp==NULL){
    fprintf(stderr,"unable to open output file \n");
    return;
  }
  fprintf(fp,"%d \n",d);
  writeNode(root,d,fp);
  fclose(fp);
}

void writeNode(treenode* node, int d, FILE* fp){
  int i;
  //n, mu, mup, R, isLeaf 
  fprintf(fp,"%d %f %d\n",node->n,node->R,node->isLeaf);  
  for (i=0;i<d;i++)
    fprintf(fp,"%f ",node->mu[i]);
  fprintf(fp,"\n");
  for (i=0;i<d;i++)
    fprintf(fp,"%f ",node->mup[i]);
  fprintf(fp,"\n");
  
  for (i=0;i<node->n;i++){
    fprintf(fp,"%d ",node->inds[i]);
  }
  fprintf(fp,"\n");
  if(!(node->isLeaf)){
    writeNode(node->l,d,fp);
    writeNode(node->r,d,fp);
  }
  
}

treenode* readTree(char*file){
  int d;
  FILE *fp = fopen(file,"r");
  if(fp==NULL){
    fprintf(stderr,"unable to open input file... exiting \n");
    exit(1);
  }
  if(fscanf(fp,"%d \n",&d)==EOF){
    fprintf(stderr, "error reading tree... exiting \n");
    exit(1);
  }
  return readNode(d,fp);
  fclose(fp);
}

treenode* readNode(int d,FILE* fp){
  int i, pos_error;
  float temp;

  treenode* node = (treenode*)calloc(1,sizeof(treenode));
  pos_error = fscanf(fp,"%d %f %d\n",&(node->n),&temp,&(node->isLeaf));
  node->R = (double)temp;
  node->mu = calloc(d,sizeof(double));
  node->mup = calloc(d,sizeof(double));
  node->inds = calloc(node->n,sizeof(int));

  for(i=0;i<d;i++){
    pos_error = fscanf(fp,"%f ",&temp);
    node->mu[i]=(double)temp;
  }
  pos_error = fscanf(fp,"\n");
  for(i=0;i<d;i++){
    pos_error = fscanf(fp,"%f ",&temp);
    node->mup[i]=(double)temp;
  }
  pos_error = fscanf(fp,"\n");
  for(i=0;i< node->n;i++){
    pos_error = fscanf(fp,"%d ",&(node->inds[i]));
  }
  pos_error = fscanf(fp,"\n");

  if( pos_error == EOF){
    fprintf(stderr, "error reading tree... exiting \n");
    exit(1);
  }
  
  if (!(node->isLeaf)){
    node->l=readNode(d,fp);
    node->r=readNode(d,fp);
  }
  return node;
}


void insert(int* NNs, double* dtoNNs, int newNN, double dNNnew, int k){
  int i=0;
  
  while (i < (k-1) && dtoNNs[i+1]>dNNnew){
    NNs[i]=NNs[i+1];
    dtoNNs[i]=dtoNNs[i+1];
    i++;
  }
  NNs[i]=newNN;
  dtoNNs[i]=dNNnew;
}

#endif

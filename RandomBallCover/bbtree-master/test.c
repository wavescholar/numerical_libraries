//This file is part of the bregman ball tree package.
//(c) copyright 2008, Lawrence Cayton

#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "bbtree.h"
#include "bregdiv.h"
#include "utils.h"
#include "search.h"

void processArgs(int,char**);
double timediff(struct timeval,struct timeval);

char *datafile, *queryfile, *outfile;
int isOutfile=0;
int n=-1,m=-1,d=-1;
double eps=0;
int divChoice=USEL2;
int bucketSize = 50;
int k=1;

int main(int argc, char** argv){
  int i,j, *NNs;
  double **x,**q, bbtime, brutetime,*dToNNs,divTemp;
  treenode *root;
  struct timeval tvB,tvE;

  printf("**** bbtree **** \n");
  processArgs(argc,argv);
  x = calloc(n,sizeof(double*));
  q = calloc(m,sizeof(double*));
  
  for(i=0;i<n;i++)
    x[i]=calloc(d,sizeof(double));
  for(i=0;i<m;i++)
    q[i]=calloc(d,sizeof(double));
  dToNNs = calloc(k,sizeof(double));
  NNs = calloc(k,sizeof(int));

  readData(x,n,d,datafile);
  readData(q,m,d,queryfile);

  bregdiv div;
  switch(divChoice){
  case USEL2:
    div = l2squared();
    printf("divergence = l_2^2\n");
    break;
  case USEKL:
    div = kl();
    printf("divergence = KL\n");
    break;
  case USEKLD:
    div = dkl();
    printf("divergence = conjugate to KL\n");
    printf("WARNING: this has not been tested thoroughly\n");
    break;
  case USEIS:
    div = is();
    printf("divergence = Itakura-Saito\n");
    printf("WARNING: this has not been tested thoroughly\n");
    break;
  }

  printf("building.....\n");
  gettimeofday(&tvB,NULL);
  root = buildTree(x,n,d,div,bucketSize);
  gettimeofday(&tvE,NULL);
  printf("done... build time: %6.2f \n",timediff(tvB,tvE));
  
  gettimeofday(&tvB,NULL);
  multisearch(root,q,x,div,n,d,m,eps,k,INT_MAX);
  gettimeofday(&tvE,NULL);
  bbtime = timediff(tvB,tvE);
  printf("BBTREE time elapsed = %6.3f \n",bbtime);
  
  
  //Example of how to save & retrieve a bbtree:

  /*writeTree(root,d,"treefile.txt");
  deleteTree(root);
  root = readTree("treefile.txt");

  gettimeofday(&tvB,NULL);
  multisearch(root,q,x,div,n,d,m,eps,0,1);
  gettimeofday(&tvE,NULL);
  bbtime = timediff(tvB,tvE);
  printf("BBTREE time elapsed = %6.3f \n",bbtime);
  */

  
  //brute force
  double curmin;
  int curBest;
  
  gettimeofday(&tvB,NULL);
  for(i=0;i<m;i++){
    for(j=0;j<k;j++){
      dToNNs[j]=HUGE_VAL;
      NNs[j]=-1;
    }
    curmin=HUGE_VAL;
    curBest=-1;
    for(j=0;j<n;j++){
      divTemp = div.div(x[j],q[i],d);
      if(NNs[0]==-1 || divTemp < dToNNs[0]){
	insert(NNs,dToNNs,j,divTemp,k);  
      }
    }
    /*    printf("query %d nns are \n",i);
    for(j=0;j<k;j++)
      printf("%d ",NNs[j]);
      printf("\n"); */
  }
  gettimeofday(&tvE,NULL);
  brutetime = timediff(tvB,tvE);
  printf("BRUTE time elapsed = %6.3f \n",brutetime);
  
  if(isOutfile){
    writeDoubs(2,outfile,bbtime,brutetime);
  }

  for(i=0;i<n;i++)
    free(x[i]);
  for(i=0;i<m;i++)
    free(q[i]);
  free(x);
  free(q);
  free(NNs);
  free(dToNNs);
  deleteTree(root);

  return 0;
}

void processArgs(int argc, char**argv){
  int i=1;
  if(argc <= 1){
    printf("usage:\n testBBT -f dataFile -q queryFile -n numPts -m numQueries -d dim [-e epsilon] [-b bregdiv] [-o outFile] [-s bucketsize] [-k numNeighbors] \n");
    printf("\n\t where bregdiv = %d for l_2^2 (default) \n\t       bregdiv = %d for KL \n\t       bregdiv = %d for exponential (conj to KL) \n\t       bregdiv = %d for Itakura-Saito\n\n", USEL2, USEKL,USEKLD,USEIS);
    exit(0);
  }
  
  while(i<argc){
    if(!strcmp(argv[i], "-d"))
      d = atoi(argv[++i]);
    else if(!strcmp(argv[i], "-n"))
      n=atoi(argv[++i]);
    else if(!strcmp(argv[i], "-m"))
      m=atoi(argv[++i]);
    else if(!strcmp(argv[i], "-f"))
      datafile = argv[++i];
    else if(!strcmp(argv[i], "-q"))
      queryfile= argv[++i];
    else if(!strcmp(argv[i], "-s"))
      bucketSize=atoi(argv[++i]);
    else if(!strcmp(argv[i], "-o")){
      outfile= argv[++i];
      isOutfile=1;
    }
    else if(!strcmp(argv[i], "-e"))
      eps=atof(argv[++i]);
    else if(!strcmp(argv[i], "-b"))
      divChoice=atoi(argv[++i]);
    else if(!strcmp(argv[i], "-k"))
      k=atoi(argv[++i]);
    else{
      fprintf(stderr,"unrecognized option.. exiting \n");
      exit(1);
    }
    i++;
  }
  
  if(n==-1 || m==-1 || d==-1 || datafile==NULL || queryfile==NULL){
    fprintf(stderr,"more arguments needed.. exiting \n");
    exit(1);
  }

  if(divChoice!=USEL2 && divChoice!=USEKL){
    fprintf(stderr,"not a legit bregman divergence.. exiting \n");
    exit(1);
  }

  if(eps<0){
    fprintf(stderr,"epsilon must be >= 0.. exiting \n");
    exit(1);
  }

}  

double timediff(struct timeval start, struct timeval end){
  return (double)(end.tv_sec+end.tv_usec/1e6 - start.tv_sec - start.tv_usec/1e6); 
}

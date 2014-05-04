/* This file is part of the Random Ball Cover (RBC) library.
 * (C) Copyright 2011, Lawrence Cayton [lcayton@tuebingen.mpg.de]
 */
#include<omp.h>
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>
#include<math.h>
#include "defs.h"
#include "utils.h" 
#include "brute.h"
#include "rbc.h"

void parseInput(int,char**);
void readData(char*,matrix);
void readDataText(char*,matrix);
void writeNeighbs(char*,char*,unint**,real**);

char *dataFileX, *dataFileQ, *dataFileXtxt, *dataFileQtxt, *outFile, *outFiletxt;
unint n=0, m=0, d=0, numReps=0, runBrute=0;
unint K=1;

int main(int argc, char**argv){
  matrix x, q;
  unint i;
  struct timeval tvB,tvE;

  printf("********************************\n");
  printf("RANDOM BALL COVER -- CPU version\n");
  printf("********************************\n");
  
  parseInput(argc,argv);

  int threadMax = omp_get_max_threads();
  printf("number of threads = %d \n",threadMax);
  
  initMat( &x, n, d );
  initMat( &q, m, d );
  x.mat = (real*)calloc( sizeOfMat(x), sizeof(*(x.mat)) );
  q.mat = (real*)calloc( sizeOfMat(q), sizeof(*(q.mat)) );
  if( !x.mat || !q.mat ){
    fprintf(stderr, "memory allocation failure .. exiting \n");
    exit(1);
  }

  
  // Read data:
  if(dataFileXtxt)
    readDataText(dataFileXtxt, x);
  else
    readData(dataFileX, x);
  if(dataFileQtxt)
    readDataText(dataFileQtxt, q);
  else
    readData(dataFileQ, q);
 
  
  unint **NNs = (unint**)calloc( m, sizeof(*NNs) ); //indices of NNs
  real **dNNs = (real**)calloc( m, sizeof(*dNNs) ); //dists to NNs
  for(i=0;i<m; i++){
    NNs[i] = (unint*)calloc(K, sizeof(**NNs));
    dNNs[i] = (real*)calloc(K, sizeof(**dNNs) );
  }

  matrix rE; //matrix of representatives
  rep *riE = (rep*)calloc( CPAD(numReps), sizeof(*riE) ); //data struct for RBC
  
  // ******** builds the RBC
  gettimeofday(&tvB,NULL);
  buildExact(x, &rE, riE, numReps);
  gettimeofday(&tvE,NULL);
  double buildTime =  timeDiff(tvB,tvE);
  printf("exact build time elapsed = %6.4f \n", buildTime );
  

  // ******** queries the RBC
  gettimeofday(&tvB,NULL);
  if( threadMax<8 )
    searchExactK(q, x, rE, riE, NNs, dNNs, K);
  else
    searchExactManyCoresK(q, x, rE, riE, NNs, dNNs, K);
  gettimeofday(&tvE,NULL);
  double searchTime =  timeDiff(tvB,tvE);
  printf("exact k-nn search time elapsed = %6.4f \n", searchTime );


  // ******** runs brute force search.
  if(runBrute){
    printf("running brute force search..\n");
    unint **NNsBrute = (unint**)calloc( m, sizeof(*NNsBrute) );
    real **dToNNsBrute = (real**)calloc( m, sizeof(*dToNNsBrute) );;
    for(i=0; i<m; i++){
      NNsBrute[i] = (unint*)calloc( K, sizeof(**NNsBrute) );
      dToNNsBrute[i] = (real*)calloc( K, sizeof(**dToNNsBrute) );
    }
    gettimeofday(&tvB,NULL);
    bruteK(x,q,NNsBrute,dToNNsBrute,K);
    gettimeofday(&tvE,NULL);
    double bruteTime = timeDiff(tvB,tvE);
    printf("brute time elapsed = %6.4f \n", bruteTime );
    
    free(NNsBrute);
    free(dToNNsBrute);
  }
  
  if( outFile || outFiletxt )
    writeNeighbs( outFile, outFiletxt, NNs, dNNs );
  
  //clean up
  freeRBC( rE, riE );

  for(i=0;i<m; i++){
    free(NNs[i]);    free(dNNs[i]);
  }
  free(NNs); free(dNNs);
  free(x.mat);
  free(q.mat);
  
  return 0;
}


void parseInput(int argc, char **argv){
  int i=1;
  if(argc <= 1){
    printf("\nusage: \n  exactRBC -x datafileX -q datafileQ -n numPts (DB) -m numQueries -d dim -r numReps [-o outFile] [-b] [-k neighbs]\n\n");
    printf("\tdatafileX    = binary file containing the database\n");
    printf("\tdatafileQ    = binary file containing the queries\n");
    printf("\tnumPts       = size of database\n");
    printf("\tnumQueries   = number of queries\n");
    printf("\tdim          = dimensionality\n");
    printf("\tnumReps      = number of representatives\n");
    printf("\toutFile      = binary output file (optional)\n");
    printf("\tneighbs      = num neighbors (optional; default is 1)\n"); 
    printf("\n\tuse -b option to run brute force search (in addition to the RBC)\n");
    printf("\n\n\tTo input/output data in text format (instead of bin), use the \n\t-X and -Q and -O switches in place of -x and -q and -o (respectively).\n");
    printf("\n\n");
    exit(0);
  }
  
  while(i<argc){
    if(!strcmp(argv[i], "-x"))
      dataFileX = argv[++i];
    else if(!strcmp(argv[i], "-q"))
      dataFileQ = argv[++i];
    else if(!strcmp(argv[i], "-X"))
      dataFileXtxt = argv[++i];
    else if(!strcmp(argv[i], "-Q"))
      dataFileQtxt = argv[++i];
    else if(!strcmp(argv[i], "-n"))
      n = atoi(argv[++i]);
    else if(!strcmp(argv[i], "-m"))
      m = atoi(argv[++i]);
    else if(!strcmp(argv[i], "-d"))
      d = atoi(argv[++i]);
    else if(!strcmp(argv[i], "-r"))
      numReps = atoi(argv[++i]);
    else if(!strcmp(argv[i], "-b"))
      runBrute = 1;
    else if(!strcmp(argv[i], "-o"))
      outFile = argv[++i];
    else if(!strcmp(argv[i], "-O"))
      outFiletxt = argv[++i];
    else if(!strcmp(argv[i], "-k"))
      K = atoi(argv[++i]);
    else{
      fprintf(stderr,"%s : unrecognized option.. exiting\n",argv[i]);
      exit(1);
    }
    i++;
  }

  if( !n || !m || !d || !numReps ){
    fprintf(stderr,"more arguments needed.. exiting\n");
    exit(1);
  }
  if( (!dataFileX && !dataFileXtxt) || (!dataFileQ && !dataFileQtxt) ){
    fprintf(stderr,"more arguments needed.. exiting\n");
    exit(1);
  }
  if( (dataFileX && dataFileXtxt) || (dataFileQ && dataFileQtxt) ){
    fprintf(stderr,"you can only give one database file and one query file.. exiting\n");
    exit(1); 
  }
  if(numReps>n){ 
    fprintf(stderr,"can't have more representatives than points.. exiting\n"); 
    exit(1); 
  } 
  if(K>n){
    fprintf(stderr,"K can't be greater than the number of points.. exiting\n"); 
    exit(1); 
  }
}


// loads data from a binary file into data in row-major order.
void readData(char *dataFile, matrix x){
  FILE *fp;
  unint numRead;
  unint i;

  fp = fopen(dataFile,"r");
  if(fp==NULL){
    fprintf(stderr,"error opening file.. exiting\n");
    exit(1);
  }
    
  for( i=0; i<x.r; i++ ){ //can't load everything in one fread
                           //because matrix is padded.
    numRead = fread( &x.mat[IDX(i, 0, x.ld)], sizeof(real), x.c, fp);
    if(numRead != x.c){
      fprintf(stderr,"error reading file.. exiting \n");
      exit(1);
    }
  }
  fclose(fp);
}


void readDataText(char *dataFile, matrix x){
  FILE *fp;
  real t;
  int i,j;

  fp = fopen(dataFile,"r");
  if(fp==NULL){
    fprintf(stderr,"error opening file.. exiting\n");
    exit(1);
  }
    
  for(i=0; i<x.r; i++){
    for(j=0; j<x.c; j++){
      if(fscanf(fp,"%f ", &t)==EOF){
	fprintf(stderr,"error reading file.. exiting \n");
	exit(1);
      }
      x.mat[IDX( i, j, x.ld )]=(real)t;
    }
  }
  fclose(fp);
}


void writeNeighbs(char *file, char *filetxt, unint **NNs, real **dNNs){
  unint i,j;
  
  if( filetxt ) { //write text

    FILE *fp = fopen(filetxt,"w");
    if( !fp ){
      fprintf(stderr, "can't open output file\n");
      return;
    }
    
    for( i=0; i<m; i++ ){
      for( j=0; j<K; j++ )
	fprintf( fp, "%u ", NNs[i][j] );
      fprintf(fp, "\n");
    }
    
    for( i=0; i<m; i++ ){
      for( j=0; j<K; j++ )
	fprintf( fp, "%f ", dNNs[i][j] );
      fprintf(fp, "\n");
    }
    fclose(fp);
    
  }

  if( file ){ //write binary

    FILE *fp = fopen(file,"wb");
    if( !fp ){
      fprintf(stderr, "can't open output file\n");
      return;
    }
    
    for( i=0; i<m; i++ )
      fwrite( NNs[i], sizeof(*NNs[i]), K, fp );
    for( i=0; i<m; i++ )
      fwrite( dNNs[i], sizeof(*dNNs[i]), K, fp );
    
    fclose(fp);
  }
}

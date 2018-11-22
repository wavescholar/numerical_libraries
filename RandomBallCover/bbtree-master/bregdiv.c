#ifndef BREGDIV_C
#define BREGDIV_C

//This file is part of the bregman ball tree package.
//(c) copyright 2008, Lawrence Cayton

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include "bregdiv.h"
#include "utils.h"

/* **************************************** */
/* GENERALIZED KL-DIVERGENCE */
/* **************************************** */
bregdiv kl(){
  bregdiv kldiv;
  kldiv.div = calckldiv;
  kldiv.gradf = calcklgradf;
  kldiv.gradfstar = calcklgradfstar;
  kldiv.divmat = calckldivmat;
  return kldiv;
}

double calckldiv(double *x, double *y, int d){
  int i;
  double ans=0;
  for(i=0;i<d;i++){
    if(x[i] > ERRORTOL )
      ans+= x[i]*log(x[i]/y[i])+y[i]-x[i];
  }
  return ans;
}

void calcklgradf(double *gf, double *x, int d){
  int i;
  for(i=0;i<d;i++)
    gf[i]=log(x[i]) +1.0;
} 

void calcklgradfstar(double *gfs,double *x, int d){
  int i;
  for(i=0;i<d;i++)
    gfs[i] = exp(x[i]-1.0);
}

void calckldivmat(double **dmat, double **x, double** y, int n, int m, int d){
  int i,j;
  for (i=0;i<n;i++){
    for(j=0;j<m;j++)
      dmat[i][j] = calckldiv(x[i],y[j],d);
  }
}


/* **************************************** */
/* l_2^2 */
/* **************************************** */
bregdiv l2squared(){
  bregdiv l2;
  l2.div = calcl2div;
  l2.gradf = calcl2gradf;
  l2.gradfstar = calcl2gradfstar;
  l2.divmat = calcl2divmat;
  return l2;
}

double calcl2div(double* x, double* y, int d){
  int i;
  double ans=0;
  for(i=0;i<d;i++)
    ans+=(x[i]-y[i])*(x[i]-y[i]);
  return ans*.5;
}

  
void calcl2gradf(double* gf, double *x, int d){
  int i;
  for(i=0;i<d;i++)
    gf[i]=x[i];
}

void calcl2gradfstar(double* gfs, double *x, int d){
  int i;
  for(i=0;i<d;i++)
    gfs[i]=x[i];
}

void calcl2divmat(double** dmat, double** x, double** y, int n, int m, int d){
  int i,j,k;

  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      dmat[i][j]=0;
      for(k=0;k<d;k++){
	dmat[i][j]+=.5*(x[i][k] - y[j][k])*(x[i][k] - y[j][k]);
      }
    }
  }

}
  

/* **************************************** */
/* dual KL */
/* **************************************** */
bregdiv dkl(){
  bregdiv dkl;
  dkl.div = calcdkldiv;
  dkl.gradf = calcdklgradf;
  dkl.gradfstar = calcdklgradfstar;
  dkl.divmat = calcdkldivmat;
  return dkl;
}



double calcdkldiv(double *x, double *y, int d){
  int i;
  double ans=0;
  for(i=0;i<d;i++)
    ans+=exp(x[i]-1) - (1+x[i]-y[i])*exp(y[i]-1);
  return ans;
}

void calcdklgradf(double *gf, double*x, int d){
  int i;
  for(i=0;i<d;i++)
    gf[i]=exp(x[i]-1);

}

void calcdklgradfstar(double* gfs, double *x, int d){
  int i;
  for (i=0;i<d;i++)
    gfs[i]=log(x[i]) +1;
}

void calcdkldivmat(double** dmat,double** x, double **y, int n, int m, int d){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      dmat[i][j]=calcdkldiv(x[i],y[j],d);
    }
  }
}

/* **************************************** */
/* Itakura-Saito */
/* **************************************** */

bregdiv is(){
  bregdiv is;
  is.div = calcisdiv;
  is.gradf = calcisgradf;
  is.gradfstar = calcisgradfstar;
  is.divmat = calcisdivmat;
  return is;
}

double calcisdiv(double *x, double *y, int d){
  double ans=0;
  int i;
  for (i=0;i<d;i++)
    ans+=x[i]/y[i] - log(x[i]/y[i]) -1;
  return ans;
}


void calcisgradf(double *gf, double *x, int d){
  int i;
  for(i=0;i<d;i++)
    gf[i]=-1/x[i];
}


void calcisgradfstar(double *gfs, double *x, int d){
  int i;
  for(i=0;i<d;i++)
    gfs[i]=-1/x[i];
  
}

void calcisdivmat(double** dmat, double** x, double **y, int n, int m, int d){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      dmat[i][j]=calcisdiv(x[i],y[j],d);
    }
  }
}


#endif

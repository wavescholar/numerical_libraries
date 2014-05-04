/* This file is part of the Random Ball Cover (RBC) library.
 * (C) Copyright 2011, Lawrence Cayton [lcayton@tuebingen.mpg.de]
 */
#ifndef DEFS_H
#define DEFS_H

#include<float.h>
#include<limits.h>
#include<math.h>
#include<stdint.h>

#define FLOAT_TOL 1e-7
#define VEC_LEN 4 
// The width of the vector (sse) instructions.  Change to 2 
// if you wish to work in double precision.

#define CL 16 
// CL is for Cache Line.  This number is not actually 
// the cache line size, but is used to make the loops more cache-efficient
// in a simple way.

#define DEF_LIST_SIZE 1024

// The following macros define the distance measure.
//#define DIST(i,j) ( fabs((i)-(j)) )  // L_1
//#define DIST_EXP(x) ( (x) ) //L_1
//#define DIST_ROOT(x) ( (x) ) //L_1

// DIST returns the distance for a single coordinate.
// DIST_EXP defines the exponent applied to the summed DIST values.
// DIST_ROOT is the inverse of DIST_EXP


// L_2 versions of the above macros:
#define DIST(i,j) ( ( (i)-(j) )*( (i)-(j) ) )  
#define DIST_EXP(x) ( sqrt(x) ) 
#define DIST_ROOT(x) ( (x)*(x) ) 

// Format that the data is manipulated in:
typedef float real;
#define MAX_REAL FLT_MAX
#define MIN_REAL (-1.0*FLT_MAX)

// To switch to double precision, comment out the above 
// three lines and uncomment the following three lines.   
// You should also change the VEC_LEN to 2.

//typedef double real;
//#define MAX_REAL DBL_MAX
//#define MIN_REAL (-1.0*DBL_MAX)

#define DUMMY_IDX UINT_MAX
#define DELETED_IDX UINT_MAX-1 //used by heap

#define IDX(i,j,ld) (((i)*(ld))+(j))
// Row major indexing

#define PAD(i) ( ((i)%VEC_LEN)==0 ? (i):((i)/VEC_LEN)*VEC_LEN+VEC_LEN ) 
// increases an int to the next multiple of VEC_LEN

#define CPAD(i) ( ((i)%CL)==0 ? (i):((i)/CL)*CL+CL )
// ditto but for CL

#define DPAD(i) ( ((i)%VEC_LEN)==0 ? (i):((i)/VEC_LEN)*VEC_LEN ) 
// decreases an int to the next multiple of VEC_LEN

#define MAX(i,j) ((i) > (j) ? (i) : (j))
#define MIN(i,j) ((i) < (j) ? (i) : (j))

#define MAXI(i,j,k,l) ((i) > (j) ? (k) : (l))
// the same as max, but takes in index arguments (k,l) and returns
// the index of the max instead of the max.

#define GETBIT(i) ( 1ul<<(i) )

typedef uint32_t unint;

typedef struct {
  real *mat;
  unint r; //rows
  unint c; //cols
  unint pr; //padded rows
  unint pc; //padded cols
  unint ld; //the leading dimension (in this code, this is the same as pc)
} matrix;


typedef struct {
  char *mat;
  unint r;
  unint c;
  unint pr;
  unint pc;
  unint ld;
} charMatrix;


typedef struct {
  unint *mat;
  unint r;
  unint c;
  unint pr;
  unint pc;
  unint ld;
} intMatrix;


typedef struct { // struct for a representative point
  unint* lr; //list of owned DB points
  real* dists; //dist to the owned points
  unint len;  //length of lr
  unint start; //used in the re-ordered version of search
  real radius;
} rep;


typedef struct { //very simple list data type
  unint* x;
  unint len;
  unint maxLen;
} intList;



typedef struct { //array-based heap element
  unint id;
  real val;
} heapEl; 


typedef struct { //array-based heap
  heapEl *h;
  unint len;
} heap;
#endif

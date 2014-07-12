/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE SMatrixA.h
   Class template for the 2-dimensional discrete Laplacian on unit
   square with zero Dirichlet boundary conditions. A is the nx*nx 
   by nx*nx block tridiagonal matrix

                              | T -I          |
                              |-I  T -I       |
                          A = |   -I  T       |
                              |        ...  -I|
                              |           -I T|

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef SMATRIXA_H
#define SMATRIXA_H

#include "matprod.h"

template<class T>
class SymMatrixA: public MatrixWithProduct<T> {

 private:

  int nx;

  void MultTv(T* x, T* y);

 public:

  void MultMv(T* v, T* w);

  SymMatrixA(int nxval): MatrixWithProduct<T>(nxval*nxval) { nx = nxval; }

}; // SymMatrixA.


template<class T>
void SymMatrixA<T>::MultTv(T* x, T* y)
/*
  Computes the matrix vector multiplication y<---T*x
  where T is a nx by nx tridiagonal matrix with DD on the
  diagonal, DL on the subdiagonal, and DU on the superdiagonal.
*/

{

  int  j;
  T    dd, dl, du;

  const T one  = 1.0;
  const T four = 4.0;

  dd  = four;
  dl  = -one;
  du  = -one;

  y[0] =  dd*x[0] + du*x[1];
  for (j = 1; j <= nx-2; j++) {
    y[j] = dl*x[j-1] + dd*x[j] + du*x[j+1];
  }
  y[nx-1] =  dl*x[nx-2] + dd*x[nx-1];

  return;

} // MultTv


template<class T>
void SymMatrixA<T>::MultMv(T* v, T* w)
/*
  Computes w <- A*v.
*/

{

  int  j, lo, n2;
  T    h2;

  const T one = 1.0;

  MultTv(&v[0],&w[0]);
  axpy(nx, -one, &v[nx], 1, &w[0], 1);

  for (j = 2; j<=nx-1; j++) {
    lo = (j-1)*nx;
    MultTv(&v[lo], &w[lo]);
    axpy(nx, -one, &v[lo-nx], 1, &w[lo], 1);
    axpy(nx, -one, &v[lo+nx], 1, &w[lo], 1);
  }

  lo = (nx-1)*nx;
  MultTv(&v[lo], &w[lo]);
  axpy(nx, -one, &v[lo-nx], 1, &w[lo], 1);

  // Scaling the vector w by (1/h^2), where h is the mesh size.

  n2 = nx*nx;
  h2 = T((nx+1)*(nx+1));
  scal(n2, h2, w, 1L);

  return;

} //  MultMv.


#endif // SMATRIXA_H

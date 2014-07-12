/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE SMatrixC.h
   Class template for the one dimensional discrete Laplacian on
   the interval [0,1] with zero Dirichlet boundary conditions.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef SMATRIXC_H
#define SMATRIXC_H

#include "matprod.h"
#include "blas1c.h"
#include "lapackc.h"

template<class T>
class SymMatrixC: public MatrixWithProduct<T> {

 public:

  void MultMv(T* v, T* w);

  SymMatrixC(int nv): MatrixWithProduct<T>(nv) { }

  virtual ~SymMatrixC() { }

}; // SymMatrixC.


template<class T>
void SymMatrixC<T>::MultMv(T* v, T* w)
/*
  Matrix-vector product.
*/

{

  int  j;
  T    h2;

  const T two = 2.0;

  w[0] =  two*v[0] - v[1];
  for (j=1; j<ncols()-1; j++) {
    w[j] = - v[j-1] + two*v[j] - v[j+1];
  }
  w[ncols()-1] = - v[ncols()-2] + two*v[ncols()-1];

  // Scaling the vector w by (1 / h).

  h2 = T(ncols()+1);
  scal(ncols(), h2, w, 1L);

  return;

} //  MultMv.


#endif // SMATRIXC_H


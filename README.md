Packages
========

The Packages git repository consists of the 
static numerical libraries used in other kl 
projects.

We maintain our own copies for these
since the port to the Intel compiler and changes to
the C++ language specification required modification 
of some source files.

Expokit

SDPA

SuiteSparse

Super LU

FastGauss

MatLib

BullMoutain Entropy RNG

Random Ball Cover

The solution for these resides in the kl root folder which is not
in source control.  The individual project files are in the repo.

The directory structure of Arpack libraries is suboptimal. ARPACK2 contains the Visual Studio Project files, arpack++\ARPACK contains the Fortran code. The C++ interface for Arpack is in arpack++\include.  Note that Arpack++ uses SuperLU and UMFPACK in the sparse examples. 

\section{ARPACK}
ARPACK++ is an object-oriented version of the Fortran ARPACK package. ARPACK is designed to compute a few eigenvalues and eigenvectors of large scale sparse matrices and pencils via the Arnoldi process for finding eigenvalues called. These methods utilize Krylov Subspace Projections for iterative solution that avoids matrix multiplication.  ARPACK implements the implicit restarted Arnoldi method which reduces the storage requirements of the traditional Lanczos iteration for Hermitian matrices and Arnoldi iteration for general matrices.  The key to the Krylov method is to calculate the linear subspace of $\Real^{(n,n)}$ induced by span of the first m powers of the image of $b$ under a linear operator $A$, $\kappa_m(A,b) | A \in \mathbb R^{(n,n)}
b\ in \mathbb R^n = \{b, Ab (A)^2b, \ldots (A)^mb \}$.  This avoids direct matrix matrix operations when finding the first few eigenvector, eigenvalue pairs in a large system of linear
equations.

\section{ATLAS}
Automatically Tuned Linear Algebra software.

\section{METIS}
METIS is a software library for finite element analysis and graph partitions.  It also can be used to reduce the fill order of
sparse matrices.


\section{SDPA}
SDPA is a software library for solving SDPs using on the Mehrotra-type predictor-corrector infeasible primal-dual interior-point method. It is implemented C++ language and utilizes the machine dependent BLAS such as Intel MKL, ATLAS. LAPACK routines are used for matrix computations.  Efficient methods to compute the search directions exploiting the sparsity of the data matrices are implemented. Sparse or dense Cholesky factorization for the Schur complemetn matrix is automatically selected. The calculation of the Schur complement
matrix is implemented in reentrant code. A sparse version of SDPA is available that uses METIS and SPOOLES libraries for finding a proper sparse structure of the problem.

\section{SPOOLS}
SPOOLES is a library for solving sparse real and complex linear systems of equations. SPOOLES can factor and solve square linear systems of equations with symmetric structure, and it can compute multiple minimum degree, generalized nested dissection and multisection orderings of matrices with symmetric structure.  SPOOLES utilizes a variety of Krylov iterative methods. The preconditioner is a drop tolerance factorization.

\section{SuperLU}
SuperLU ( http://crd-legacy.lbl.gov/~xiaoye/SuperLU/) is a general purpose library for the direct solution of large, sparse, nonsymmetric systems of linear equations on high performance machines. The library is written in C and is callable from either C or Fortran. The library routines will perform an LU decomposition with partial pivoting and triangular system solves through forward and back substitution. The LU factorization routines can handle non-square matrices but the triangular solves are performed only for square matrices. The matrix columns may be preordered (before factorization) either through library or user supplied routines. This pre-ordering for sparsity is completely separate from the factorization. Working precision iterative refinement subroutines are provided for improved backward stability. Routines are also provided to equilibrate the system, estimate the condition number, calculate the relative backward error, and estimate error bounds for the refined solutions.

\section{SuiteSparse}
Tim Davis' ( http://www.cise.ufl.edu/~davis/welcome.html) collection of sparse matrix software.  Tim is also the curator of The University of Florida Sparse Matrix Collection (http://www.cise.ufl.edu/research/sparse/matrices/) a must see for anyone interested in sparse
matrices and visualization.

AMD: symmetric approximate minimum degree
BTF: permutation to block triangular form
CAMD: symmetric approximate minimum degree
CCOLAMD: constrained column approximate minimum degree
COLAMD: column approximate minimum degree
CHOLMOD: sparse supernodal Cholesky factorization and update/downdate
CSparse: a concise sparse matrix package
CXSparse: an extended version of CSparse
KLU: sparse$ LU$ factorization, for circuit simulation
LDL: a simple $LDL^T$ factorization
UMFPACK: sparse multifrontal $LU$ factorization
RBio: MATLAB toolbox for reading/writing sparse matrices
UFconfig: common configuration for all but CSparse
SuiteSparseQR: multifrontal sparse $QR$

\subsection{AMD}
AMD is a set of routines for pre-ordering a sparse matrix prior to numerical factorization. It uses an approximate minimum degree ordering algorithm to find a permutation matrix P so that the Cholesky factorization $PAP^\dag =LL^\dag$ has fewer (often much fewer) nonzero entries than the Cholesky factorization of A. The algorithm is typically much faster than other ordering methods and minimum degree ordering algorithms that compute an exact degree . Some methods, such as approximate deficiency [Rothberg and Eisenstat 1998] and graph-partitioning based methods [Hendrickson and Rothberg 1999; Karypis and Kumar 1998; Pellegrini et al. 2000; Schulze 2001] can produce better orderings, depending on the matrix. The algorithm starts with an undirected graph representation of a symmetric sparse matrix . Node $i$ in the graph corresponds to row and column i of the matrix, and there is an edge $(i,j)$ in the graph if $a_{ij}$ is nonzero. The degree of a node is initialized to the number of off diagonal non-zeros in row $i$, which is the size of the set of nodes adjacent to $i$ in the graph.

\subsection{UMFPACK}
UMFPACK is a set of routines for solving systems of linear equations, $Ax = b$, when $A$ is sparse and unsymmetric. It is based on the Unsymmetric-pattern MultiFrontal method. UMFPACK factorizes $PAQ$, $PRAQ$ and $PR^{-1}AQ$, into the product $LU$, where $L$ and $U$ are lower and upper triangular, respectively, $P$ and $Q$ are permutation matrices, and $R$ is a diagonal matrix of row scaling factors (or $R = I$ if row-scaling is not used). Both $P$ and $Q$ are chosen to reduce fill-in (new nonzeros in $L$ and $U$ that are not present in $A$). The permutation $P$ has the dual role of reducing fill-in and maintaining numerical accuracy (via relaxed partial pivoting and row interchanges). The sparse matrix $A$ can be square or rectangular, singular or non-singular, and real or complex (or any combination). Only square matrices $A$ can be used to solve $Ax = b$ or related systems. Rectangular matrices can only be factorize 



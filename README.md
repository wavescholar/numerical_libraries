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

The solution for these resides in the kl root folder which is not
in source control.  The individual project files are in the repo.

The directory structure of these libraries is suboptimal.  The packages were added during an organic evolution of one of my projects a long time ago and I’ve never revisited the location of some of the libraries.  For instance Arpack is located here;  ConvexOptimization / SDPA_INTEL_BLAS / arpack++ / ARPACK.  SDPA is a convex solver.  Arpack obviously; is not. 

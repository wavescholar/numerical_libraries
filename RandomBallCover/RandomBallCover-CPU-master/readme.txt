***[cpu] Random Ball Cover (RBC) v0.1***
Lawrence Cayton
work@lcayton.com

(C) Copyright 2011, Lawrence Cayton [work@lcayton.com]
 
This program is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

---------------------------------------------------------------------
SUMMARY

This is a C implementation of the Random Ball Cover data structure for
fast nearest neighbor (NN) search, designed for shared-memory systems.
All parallelization is handled through OpenMP.  This code contains
both the one-shot (approximate) search algorithm, and the exact 
search algorithm.  

There is a different implementation available that runs on a GPU;
visit my (Lawrence Cayton's) webpage for details.

Detailed information about the theory behind this method and its
empirical performance can be found in the paper

L. Cayton. Accelerating nearest neighbor search on manycore systems. 
Twenty-Sixth IEEE International Parallel and Distributed Processing 
Symposium (IPDPS), 2012. 


---------------------------------------------------------------------
COMPILATION

This code currently requires the GNU Scientific Library (GSL), which
is available for free on the web (or through a Linux package
manager).  It also requires the OpenMP libraries and GCC.  

To build the code, type make in a shell.  

The code has been tested under Linux and MacOS X.  


---------------------------------------------------------------------
USAGE

Two sample drivers are provided, one for the exact search algorithm,
the other for the one-shot algorithm.  Type 
$ exactRBC
or
$ oneShotRBC
at the prompt to get a list of options.  

The output file format is a list of the queries' NNs,
followed by a list of the distances to those NNs.

Basic functionality is provided through these drivers, but I recommend
integrating the RBC code directly into your code for the best
results.  In particular, the best way to use the current
implementation is to build the RBC once, then query it many times.

Both methods require the setting of one parameter, the number of
representatives.  For the one-shot algorithm, this parameter allows
one to trade-off between search quality and search speed.  For the exact
search algorithm, the results will be the same regardless of this
setting, but search performance will vary.  The best way to set this
parameter is to try a few different values out; a good starting point
is generally 5*sqrt(n), where n is the number of database points.  
See the paper for detailed information on this parameter.


---------------------------------------------------------------------
FILES

* brute.{c,h} -- implementation of brute force NN search and many
  variants.  These routines perform virtually all the work.
* rbc.{c,h} -- the core RBC algorithms.  This includes build and
  search routines for exact and approximate search.  The searchExact
  method comes in two forms, one with ManyCores appended to the
  function name; see below for discussion.  
* utils.{c,h} -- supporting code, including the implementations of
  some basic data structures and various routines useful for
  debugging.  
* dists.{c,h} -- functions that compute the distance
* defs.h -- defintions of constants and macros, including the
  distance metric.


-->DRIVERS
* exactDriver.c -- example driver for the exact search algorithm
* oneShotDriver.c -- example driver for the one-shot search 
  algorithm.


---------------------------------------------------------------------
MISC NOTES ON THE CODE


* The current version requires the GNU Scientific Library (GSL).  The
  code only uses the library for sorting and random number
  generation.  This dependency will probably be removed in the future.

* For both search methods, there are separate 1-NN and K-NN functions 
  (distinguishable by the function names).  One can run the K-NN
  functions with K=1 and get the same answer as the 1-NN functions,
  but the 1-NN functions are slightly faster, and easier to
  understand.  

* There are two versions of the exact search method for the RBC,
  searchExact(..) and searchExactManyCores(..), plus K-NN versions of
  both.  searchExact(..) is somewhat faster for systems with a small
  number of cores.  I recommend using searchExactManyCores(..) for
  systems with more than 4 cores, though you might try both methods.

* The algorithms implemented here work for an arbitrary metric.  The
  implementation currently only supports the L_1 and L_2 distance.  An
  arbitrary L_p distance is trivial to add by redefining 3 constants.
  See defs.h for an example.  If you wish to implement your own
  metric, simply replace the two functions in dists.{c,h}.  

* The code uses the default number of threads defined by OpenMP.
  Generally, this will be the number of cores in your computer.  If
  you wish to manually set this, set the OMP_NUM_THREADS environment
  variable before calling the driver.  

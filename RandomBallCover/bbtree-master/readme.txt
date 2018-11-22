***Bregman Ball Trees (bbtrees)***
Lawrence Cayton
work@lcayton.com

(C) Copyright 2008, Lawrence Cayton

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
with this program.  If not, see <http://www.gnu.org/licenses/>.

-----------------------------------------------------------

This is a C implementation of the bregman ball tree data structure 
described in 

L. Cayton, Fast nearest neighbor retrieval for bregman divergences.  
Twenty-Fifth International Conference on Machine Learning (ICML), 2008.

The code provides an implementation of the build, search, and 
approximate search algorithms described in the paper.   The software
currently supports the following divergences:
	  * l_2^2
	  * KL-divergence (relative entropy)
	  * The conjugate to the KL-divergence (exponential)
	  * Itakura-Saito divergence
It is simple to extend to other divergences.

The code is not seriously optimized; it can probably be sped up
with a bit of programming expertise.  

-----------------------------------------------------------
FILES

* bbtree.{h,c}  -- contains the build algorithm and the definition of
the bbtree data structure.

* bregdiv.{h,c} -- defines the bregdiv data structure and provides
implementations of the supported divergences.  The bregdiv data
structure is simply a struct with function pointers to the appropriate
divergence functions.

* search.{h,c} -- contains the search algorithms.

* test.c -- code for experimenting with the bbtree data
structure.  

* utils.{h,c} -- contains procedures for saving and loading bbtrees,
and a variety of other utilities used by the bbtree software.  All
constants are defined here.

-----------------------------------------------------------
COMPILING

Type make in a shell.  Compilation requires GCC.

-----------------------------------------------------------
USING 

To deploy the bbtree data structure in an application, you will likely
need to integrate this code into yours.  The test.c file provides an
example of how to use the software.  To experiment with it, type 
$ testBBT
at the prompt and a list of options will be displayed.  

The default bucketsize is 50, which was used in the exact NN
experiments in the paper.  A bucketsize of 10 was used for the 
approximate NN experiments in the paper---this is slower and was
only used so that the tradeoff between execution time and solution
quality could be examined closely.  

Approximate NN search can be accomplished in two ways. First, one can 
pass epsilon > 0 to the search procedure.  The returned point is 
guaranteed to be a (1+epsilon)-approximate NN---i.e. a point x is 
returned that satisfies
   d(x,q) <= (1+epsilon)d(x_q,q),
where x_q is the true NN.  This functionality was not explored in the 
paper.  The second method is by early stopping.  To use this feature, 
one must specify a bound on the number of leaves that will be
explored.  The approximate NN plots in the paper were produced by 
varying the leaf bound from 2^0 to 2^{log n} (i.e., full exploration)
in increments of powers of 2.  

To add your own bregman divergence, follow the prescription in
bregdiv.{h,c}.  You must create functions for computing 
* the divergence between two vectors;
* the gradient;
* the gradient of the conjugate function (ie, the inverse of the
gradient);
* the matrix of divergences between two sets of vectors.
These functions are then encapsulated in a struct of function pointers.


-----------------------------------------------------------
DATA FORMAT FOR THE TEST PROGRAM

The data must be in plain text, one vector per line, and separated by
spaces.  All points are read as doubles.  
 

-----------------------------------------------------------
ISSUES

Unbounded interpoint distances can trip up the k-means implementation
used in the build procedure.  The problem is that k-means is currently
initialized with a random point, rather than an average.


-----------------------------------------------------------
QUESTIONS, COMMENTS, ETC

Please contact work@lcayton.com with any bug reports, questions,
comments, suggestions, significant optimizations, etc.  I'm also
interested to hear about any applications of the bbtree.  

This directory contains a parallel implementation of the Parameterized
Expectations Algoritm. The code implements the algorithm described by
Maliar and Maliar in "Parameterized Expectations Algorithm and the Moving Bounds"
Journal of Business and Economic Statistics, 21, January 2003, pp. 88-92.

Matlab code by Maliar and Maliar was helpful in writing the present code, but the code
here is a new implementation with a completely different user interface, parallel
capacity, and some differences in the details of the algorithm. The Maliar and Maliar
code is available at the following sources:
http://www.amstat.org/publications/jbes/upload/index.cfm?fuseaction=ViewArticles&pub=JBES&issue=03-1-JAN
http://dge.repec.org/codes/maliar/PEA.ZIP

This Octave code parallelizes both the model simulation and the NLS
fit. The NLS fit uses analytic gradients, that are passed across the
nodes along with the objective function contributions.

To try it out on a single computer, get into Octave, and type "pea_example".

To try it out on a cluster, get into Octave and type "parallel_performance"
(this requires a running LAM, and an installation of MPITB for GNU Octave).
The ParallelKnoppix CD will help you to create a cluster. It has this code and
all dependencies installed and ready to use.

The simple_growth.oct file was compiled using Octave 2.1.73. To re-make it
with other versions, type "mkoctfile simple_growth.cc"

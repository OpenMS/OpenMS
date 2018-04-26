MIToolbox
=========
v3.0.1 for C/C++ and MATLAB/Octave

MIToolbox contains a set of functions to calculate information theoretic
quantities from data, such as the entropy and mutual information.  The toolbox
contains implementations of the most popular Shannon entropies, and also the
lesser known Renyi entropy. The toolbox also provides implementations of 
the weighted entropy and weighted mutual information from "Information Theory
with Application", S. Guiasu (1977). The toolbox only supports discrete distributions,
as opposed to continuous. All real-valued numbers will be processed by x = floor(x).

These functions are targeted for use with feature selection algorithms rather 
than communication channels and so expect all the data to be available before 
execution and sample their own probability distributions from the data.

All functions expect the inputs to be vectors or matrices of doubles.

Functions contained:
 - Entropy
 - Conditional Entropy
 - Mutual Information
 - Conditional Mutual Information
 - generating a joint variable
 - generating a probability distribution from a discrete random variable
 - Renyi's Entropy
 - Renyi's Mutual Information
 - Weighted Entropy
 - Weighted Mutual Information
 - Weighted Conditional Mutual Information

Note: all functions are calculated in log base 2, so return units of "bits".

MIToolbox works on discrete inputs, and all continuous values **must** be
discretised before use with MIToolbox. Real-valued inputs will be discretised
with x = floor(x) to ensure compatibility. MIToolbox produces unreliable
results when used with continuous inputs, runs slowly and uses much more memory
than usual. The discrete inputs should have small cardinality, MIToolbox will
treat values {1,10,100} the same way it treats {1,2,3} and the latter will be
both faster and use less memory. This limitation is due to the difficulties in
estimating information theoretic functions of continuous variables.

======

Examples:

```
>> y = [1 1 1 0 0]';
>> x = [1 0 1 1 0]';
```
```
>> mi(x,y)       %% mutual information I(X;Y)
ans =
    0.0200
```
```
>> h(x)          %% entropy H(X)
ans =
    0.9710
```
```
>> condh(x,y)    %% conditional entropy H(X|Y)
ans =
    0.9510
```
```
>> h( [x,y] )    %% joint entropy H(X,Y)
ans =
    1.9219
```
```
>> joint([x,y])  %% joint random variable XY
ans =
     1
     2
     1
     3
     4
```
======

To compile the library for use in MATLAB/OCTAVE, execute CompileMIToolbox.m
from within MATLAB in the "matlab" folder.

To compile the library for use with C programs run 'make x86' for a 32-bit
library, or 'make x64' for a 64-bit library. Then run 'sudo make install' to
install MIToolbox into /usr/local/lib & /usr/local/include.

All code is licensed under the 3-clause BSD license.

Update History
 - 08/02/2017 - v3.0.1 - Bug fix to ensure ANSI C compatibility.
 - 07/01/2017 - v3.0.0 - Refactored internals to expose integer information theoretic calculations.
 - 10/01/2016 - v2.1.2 - Relicense from LGPL to BSD. Added checks to ensure input MATLAB types are doubles.
 - 02/02/2015 - v2.1.1 - Fixed up the Makefile so it installs the headers too.
 - 22/02/2014 - v2.1  - Fixed a couple of bugs related to memory handling.
                     Added a make install for compatibility with PyFeast.
 - 30/07/2011 - v2.00 - Added implementations of the weighted entropy and weighted
                     mutual information. More cleanup of Mex entry point
                     to further check the inputs.
 - 08/11/2011 - v1.03 - Minor documentation changes to accompany the JMLR publication.
 - 15/10/2010 - v1.02 - Fixed bug where MIToolbox would cause a segmentation fault
                     if a x by 0 empty matrix was passed in. Now prints an 
                     error message and returns gracefully.
 - 02/09/2010 - v1.01 - Fixed a bug in CMIM.m where the last feature would not be 
                     selected first if it had the highest MI.
 - 07/07/2010 - v1.00 - Initial Release.
                    

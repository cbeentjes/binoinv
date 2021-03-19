# binoinv: Fast evaluation of the inverse binomial CDF

## Synopsis

The function [normcdfinv(u)](https://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__DOUBLE.html#group__CUDA__MATH__DOUBLE_1g78e93df6c3fbade8628d33e11fc94595) in NVIDIA's [CUDA maths library](https://docs.nvidia.com/cuda/cuda-math-api/index.html) (and similar functions in various libraries for CPU execution) computes the inverse of the cumulative distribution function (CDF) for the Normal distribution. It can be used to convert uniformly distributed pseudo-random (or quasi-random) numbers into psuedo-random (or quasi-random) Normals.

Here we present an article, and associated software, on algorithms for a function binoinv(u) which performs the corresponding task for binomial distributions, and can be used to generate binomial random variates.

Some aspects of the implementation, as well as the testing, build on previous work by M. B. Giles on the [inverse Poission CDF function](https://github.com/cbeentjes/poissinv).

For mathematical details, see

    @article{giles2021,
    author = {Beentjes, C. H. L., Giles, M. B.},
    title = {Algorithm XXX: approximation of the inverse binomial cumulative distribution function},
    }

For more information on fast inverse CDF computation [click here](https://people.maths.ox.ac.uk/~gilesm/codes/).

## Authors

* M. B. Giles   <mike.giles@maths.ox.ac.uk>
* C. H. L. Beentjes <casper.beentjes@maths.ox.ac.uk>

## Installation

To download:

    $ git clone https://github.com/cbeentjes/binoinv.git
    $ cd binoinv

## Matlab integration via MEX files

The fast C routines can be directly used in MATLAB applications via their MEX interfaces provided. To compile the interfaces run the MATLAB helper file compile_interface.m via the MATLAB IDE
    
    $ cd binoinv/mex
    $ compile_interface

To add CDF, inverse CDF and inverse complementary CDF MEX routines to MATLAB's search path via command line:

    $ export MATLABPATH=/path/to/binoinv/mex/:$MATLABPATH

To add said routines via MATLAB IDE:

    $ addpath(/path/to/binoinv/mex/)

To make this permanent consider [adding this line to your startup.m file](mathworks.com/help/matlab/matlab_env/add-folders-to-matlab-search-path-at-startup.html).

## Repository content

* [/analysis](https://github.com/cbeentjes/binoinv/src/master/analysis/) contains MATLAB and Mathematica code used to construct and verify parts of the codebase and code used to prepare figures for the paper.

* [/mex](https://github.com/cbeentjes/binoinv/src/master/mex/) contains MEX wrappers to call the fast CDF, inverse CDF and inverse complementary CDF routines as MATLAB routines.

* [/src](https://github.com/cbeentjes/binoinv/src/master/src/) contains several header files with the core C-routines to evaluate the inverse CDF and complementary CDF functions (binoinv.h,binocinv.h) and the regular CDF function (binocdf.h). Also present are single precision variants in the (binoinvf.h,binocinvf.h,binocdff.h).

* [/test](https://github.com/cbeentjes/binoinv/src/master/test/) contains test code to verify accuracy and test performance of the inverse CDF routines routines. The validation part uses a quad-precision function which requires gcc's [libquadmath](https://gcc.gnu.org/onlinedocs/libquadmath/) library.
 
## Licensing and acknowledgements

This code is freely available to all under a GPL license -- anyone requiring a more permissive license for commercial purposes should contact M. B. Giles.

We would be grateful if academic users would cite the paper above, once it has been published.

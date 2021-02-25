%%
% Helper file to compile the MEX interfaces for the fast inverse CDF and
% inverse complementary CDF functions for the Poisson distribution

% Double precision variants
mex binocdf_fast.c -I../src/Serial
mex binoinv_fast.c -I../src/Serial
mex binocinv_fast.c -I../src/Serial

% Single precision variants
mex binocdff_fast.c -I../src/Serial
mex binoinvf_fast.c -I../src/Serial
mex binocinvf_fast.c -I../src/Serial
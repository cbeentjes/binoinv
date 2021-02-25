% binocinvf_fast.m Help file for binocinvf_fast MEX-file.
%  binocinvf_fast.c - Fast evaluation of the Binomial inverse complementary
%  cumulative distribution function in single precision.
%
%  X=binocinvf_fast(Y,N,P) returns the inverse of the binomial 
%  complementary cdf with parameters N and P. Since the binomial 
%  distribution is discrete, binocinvf_fast returns the least integer X such
%  that the binomial complementary cdf evaluated at X, equals or exceeds Y.
%
%  The size of X is the common size of the input arguments. A scalar input  
%  functions as a constant matrix of the same size as the other inputs.    
%
%  Note that X takes the values 0,1,2,...,N.
% 
%   Reference:
%      [1]  C.H.L. Beentjes and M.B. Giles, "Algorithm XXX: approximation 
%           of the inverse binomial cumulative distribution function"
%   https://github.com/cbeentjes/binoinv
% binocdf_fast.m Help file for binocdf_fast MEX-file.
%  binocdf_fast.c - Fast evaluation of the Binomial cumulative distribution
%  function.
%
%  Y=binocdf_fast(X,N,P) returns the binomial cumulative distribution
%    function with parameters N and P at the values in X.
% 
%  The size of Y is the common size of the input arguments. A scalar input  
%  functions as a constant matrix of the same size as the other inputs.    
% 
%  The algorithm uses the cumulative sums of the binomial masses in the
%  tails of the distribution and relies on the direct calculation of the
%  regularised incomplete beta function otherwise.
%  
%  Y=binocdf_fast(X,N,P,1) returns the upper tail probability of the 
%  binomial distribution with parameters N and P at the values in X.
% 
%   Reference:
%      [1]  C.H.L. Beentjes and M.B. Giles, "Algorithm XXX: approximation 
%           of the inverse binomial cumulative distribution function"
%   https://github.com/cbeentjes/binoinv
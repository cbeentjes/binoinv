%%
% Compare speed of native MATLAB implementation with using MEX version of
% binocdf.h implementation. 
% Consider both 
% 1) the case of scalar p, which should suit MATLAB approach of
%                                               vectorised summations.
% 2) the case of vector p, which causes serious speed degradation for
%        MATLAB
%
% Conclusion for 1) only for N relatively small and number of samples large
% is MATLAB vectorised approach faster. For either N > 1e5 or small number
% of samples MEX is still faster. For N > 1e5 MATLAB starts to use a
% different summation technique which needs to evaluate the binomial PDF,
% this slows down the evaluation significantly.
%
% Conclusion for 2) MEX approach is always faster.
%

addpath('../../mex/')
if ~exist('binocdf_fast')
    mex ../../mex/binocdf_fast.c -outdir ../../mex -I../../src
end
clear all
rng(1)

%% Scalar P
disp('--scalar p--')
% Number of samples
M = 1e6;

% Binomial parameters
N = 4e4;
p = 0.5;

k = randi(N,1,M);
disp('MEX binocdf:')
tic; U1 = binocdf_fast(k,N,p); toc
disp('MATLAB binocdf:')
tic; U2 = binocdf(k,N,p); toc

E = abs(U1-U2);
rE = E./U2;
ix = find(U2 > 0);
[e,jx] = max(E(ix)./U2(ix));
kx = ix(jx);
% Max relative error location and parameters
N
p
k(kx)
e

%% Mixed P
disp('--Mixed p--')
% Number of samples
M = 1e3;

% Binomial parameters
N = 1e4;
p = rand(1,M);

k = randi(N,1,M);
disp('MEX binocdf:')
tic; U1 = binocdf_fast(k,N,p); toc
disp('MATLAB binocdf:')
tic; U2 = binocdf(k,N,p); toc

E = abs(U1-U2);
rE = E./U2;
ix = find(U2 > 0);
[e,jx] = max(E(ix)./U2(ix));
kx = ix(jx);

% Max relative error location and parameters
N
p(kx)
k(kx)
e

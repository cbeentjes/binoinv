%%
% Test erfcx and erfcxf implementation accuracy for values occuring in 
% evaluation of betaexp function. Note that the erfcx argument in that case
% is always positive, sqrt(a*rlog1(-lambda/a) + b*rlog1(lambda/b). We
% actually know that the argument for erfcx lies in [0,30].
% 
% We see that the accuracy of the both single and double precision
% implementations are very accurate.
%
%%
% Add MEX versions of the C-code erfcx/erfcxf evaluation
clc
addpath('../../../mex/')
if ~exist('erfcx_fast')
cend
if ~exist('erfcxf_fast')
    mex ../../../mex/erfcxf_fast.c -outdir ../../../mex -I../../../src/float
end
%%
clear all
close all
N = 1e7;
X = linspace(0,30,N);

% Reference values in double precision (native MATLAB implementation)
yd = erfcx(X);

%% Single precision
% Matlab implementation
ys = erfcxf_implementation(single(X));
% C implementation (not vectorised so slower for large number of samples)
% ys = erfcxf_fast(single(X));

% Absolute error
figure()
subplot(2,1,1)
semilogy(X,abs(ys-yd))
hold on
semilogy(X,eps(ys))
xlabel('x')
ylabel('absolute error')
title('single precision')

% Relative error
subplot(2,1,2)
plot(X,abs(ys-yd)./eps(ys))
xlabel('x')
ylabel('relative error (ULP)')
max(abs(ys-yd)./eps(ys))

%% Double precision
% Matlab implementation
yd2 = erfcx_implementation(X);
% C implementation (not vectorised so slower for large number of samples)
% yd2 = erfcx_fast(X);

% Absolute error
figure()
subplot(2,1,1)
semilogy(X,abs(yd-yd2))
hold on
semilogy(X,eps(yd))
xlabel('x')
ylabel('absolute error')
title('double precision')

% Relative error
subplot(2,1,2)
plot(X,abs(yd2-yd)./eps(yd))
xlabel('x')
ylabel('relative error (ULP)')
max(abs(yd2-yd)./eps(yd))

%%
%   Compare the error bounds (delta) derived for the Normal asymptotic
%   approximation and the Temme uniform asymptotic approximation.
%
%   The normal asymptotic approximation error grows away from the central
%   region, i.e. U close to 0 or 1. The Temme approximation on the other
%   hand only grows slowly.
%

clear all
close all
clc
addpath('../../mex/')
if ~isfile('../../mex/binoinv_fast')
    mex ../../mex/binoinv_fast.c -outdir ../../mex -I../../src/Serial
end

p = 0.125;
N = 1e6;

W = linspace(-3,3,1e3);
U = normcdf(W);

% Delta from Normal approximation
d1 =  (W.^4/400 +    W.^2 /200  + 1/100 )*(1+p)*(2-p).* ...
    (abs(1-2*p) + 0.25./sqrt(N*p*(1-p)))./(N*p*(1-p));

% Delta from Temme approximation
xi = binoinv_fast(U,N,p)./(N+1);

d2 = 2.25e-2./((N+1)*(xi-xi.*xi));
% p-dependent improvements to delta Temme approximation
d3 = d2*(abs(1-2*p) + 20/(sqrt(N)));
d4 = d2*(abs(1-2*p)*(1- 1*20/(sqrt(N))) + 20/(sqrt(N)));

% Plot different delta forms
figure()
plot(W,d1);
hold on
plot(W,d2);
xlabel('w')
ylabel('\delta')
legend({'\delta_N','\delta_T'})
% plot(W,d3,'o');
% plot(W,d4,'.');

% Ratio of Normal delta to Temme delta
figure()
plot(W,d1./d2)
xlabel('w')
ylabel('\delta_N / \delta_T')
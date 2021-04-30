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

lg1 = {};
lg2 = {};

%% Binomial distribution
addpath('../../mex/')
if ~isfile('../../mex/binoinv_fast.mexa64')
    mex ../../mex/binoinv_fast.c -outdir ../../mex -I../../src/Serial
end

p = 0.125;
N = 1e4;

W = 3*linspace(-1,1,1e3);
U = normcdf(W);

% Delta from Normal approximation
d1 =  (W.^4/400 +    W.^2 /200  + 1/100 )*(1+p)*(2-p).* ...
    (abs(1-2*p) + 0.25./sqrt(N*p*(1-p)))./(N*p*(1-p));

% Delta from Temme approximation
xi = binoinv_fast(U,N,p)./(N+1);

d2 = 2.25e-2./((N+1)*(xi-xi.*xi));

% Make sure to restrict to only values of W for which the asymptotic is
% actually used
ix = and((N+1)*xi >= 9, (N+1)*xi <= N - 9);
W = W(ix);
d1 = d1(ix);
d2 = d2(ix);

% Plot different delta forms
figure(1)
plot(W,d1);
hold on
plot(W,d2);
xlabel('w')
ylabel('\delta')
lg1 = [lg1(:)', {'\delta_N (binomial)'},{'\delta_T (binomial)'}];
legend(lg1)

% Ratio of Normal delta to Temme delta
figure(2)
plot(W,d1./d2)
hold on
xlabel('w')
ylabel('\delta_N / \delta_T')
lg2 = [lg2(:)', {'\delta_N / \delta_T (binomial)'}];
legend(lg2)

% Check the interval in which the ratio falls for |w| < 3
[min(d1./d2),max(d1./d2)]

%% Poisson distribution
% lam = 100;
lam = N*p; % same mean as binomial distribution for comparison

W = 3*linspace(-1,1,1e3);
U = normcdf(W);

% Delta from Normal approximation
d1 =  (W.^4/160 +    W.^2 /80  + 1/40 )./lam;

% Delta from Temme approximation
r = poissinv(U,lam)./lam;

d2 = 0.01./(lam.*r);

% Make sure to restrict to only values of W for which the asymptotic is
% actually used
ix = lam.*r >= 9;
W = W(ix);
d1 = d1(ix);
d2 = d2(ix);

% Plot different delta forms
figure(1)
plot(W,d1);
hold on
plot(W,d2);
xlabel('w')
ylabel('\delta')
lg1 = [lg1(:)', {'\delta_N (Poisson)'},{'\delta_T (Poisson)'}];
legend(lg1)

% Ratio of Normal delta to Temme delta
figure(2)
plot(W,d1./d2)
xlabel('w')
ylabel('\delta_N / \delta_T')
lg2 = [lg2(:)', {'\delta_N / \delta_T (Poisson)'}];
legend(lg2)


% Check the interval in which the ratio falls for |w| < 3
[min(d1./d2),max(d1./d2)]

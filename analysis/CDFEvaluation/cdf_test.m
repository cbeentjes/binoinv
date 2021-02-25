%%
%  Evaluate accuracy of the Binomial CDF evaluation procedures
%  Compare both against MATLAB's binocdf (direct summation)
%  and possibly against Mathematica arbitrary precision output. To get
%  Mathematica output run BinoCDFeval.nb notebook with same parameters as
%  this script
%  
%  Note that the output of binom_cdf or binocdf_fast is more accurate than
%  the binocdf code in the tails of the distribution due to inaccuracies in
%  binopdf native MATLAB code.
%
%%
% Add MEX version of the C-code CDF evaluation
clc
addpath('ERFCX/')
addpath('../../mex/')
if ~isfile('../../mex/binocdf_fast.mexa64')
    mex ../../mex/binocdf_fast.c -outdir ../../mex -I../../src/Serial
end
%%
clear all
close all

N = 2^20;
N = 52;
% Make sure to choose p that does not suffer from floating point 
% representation error! Must take p = i/2^j for 0<=i<=2^j, i,j pos. int.
i = 1;
j = 1;
p = i/2^j
q = 1-p;

% nn = 100:2100;
% nn = 0:N;
nn = floor(max(N*p - 37.5*sqrt(N*p*(1-p)),0)):ceil(min(N,N*p + 37.5*sqrt(N*p*(1-p))));

% CDF value array
Ch = zeros(2,length(nn));
CM = zeros(1,length(nn));
% PMF value array
Ph = zeros(1,length(nn));
PM = zeros(1,length(nn));

% Matlab implementation of the C-code CDF and PMF evaluation
tic
for ix = 1:length(nn)
    k = nn(ix);
    [a,b] = binom_cdf(k,N,p,q);    
    if k > N*p
        Ch(1,ix) = b;
    else
        Ch(1,ix) = a;
    end
    Ph(1,ix) = binom_pmf(k,N,p,1-p);
end
toc
% Matlab native Binomial CDF evaluation
CM(1,nn > N*p)  = binocdf(nn(nn >  N*p),N,p,'upper');
CM(1,nn <= N*p) = binocdf(nn(nn <= N*p),N,p);

% MEX version of C-code CDF evaluation
Ch(2,nn > N*p)  = binocdf_fast(nn(nn >  N*p), N, p, 1);
Ch(2,nn <= N*p) = binocdf_fast(nn(nn <= N*p), N, p);

% Matlab native PMF evaluation
PM(1,:) = binopdf(nn,N,p);

% Relative error to MATLAB evaluation
errM = Ch - CM;
errM = errM./CM;

% Relative error binom_cdf MATLAB vs MEX
errC = Ch(1,:) - Ch(2,:);
errC = errC./Ch(1,:);

pltx = 1;
figure(1)
hold on

plt(pltx) = plot(nn,abs(errM(1,:)),'.');
lg{pltx} = 'Relative error binomcdf and binocdf';
pltx = pltx + 1;

% plt(pltx) = plot(nn,abs(errM(2,:)));
% lg{pltx} = 'Relative error C-code and binocdf';
% pltx = pltx + 1;

% Plot as reference line the relative machine precision of the
% log-probability. Note that most routines calculate the CDF or PMF first
% on the log-scale so the accuracy is limited by the machine precision of
% this value.
U = Ch(1,:);
plt(pltx) = plot(nn,eps*abs(log(U)),'k','LineWidth',2);
lg{pltx} =  'eps * log-CDF';
pltx = pltx + 1;

plt(pltx) = plot(nn,eps(abs(log(U))),'LineWidth',2);
lg{pltx} = 'Machine precision log CDF';
pltx = pltx + 1;

% Boundary at which binom_cdf / MEX code starts to use direct summation
plot((N+1)*p/(2-p)*[1,1],[min(abs(errM(abs(errM)>0))),1],'r')
plot((-1+(N+1)*2*p/(1+p))*[1,1],[min(abs(errM(abs(errM)>0))),1],'r')
set(gca,'YScale','log')
legend(plt, lg(:))

%% Load more accurate data from Mathematica
% Take as reference values the arbitrary precision values from Mathematica
load('qdata.mat');
c1 = C(:,1)';
c2 = C(:,2)';
nn = C(:,3)';
c = c1;
ix = nn > N*p;
c(ix) = c2(ix);

err1 = Ch(1,:) - c; % Error MATLAB implementation binom_cdf
err2 = Ch(2,:) - c; % Error C-code implementation binom_cdf
err3 = CM(1,:) - c; % Error native MATLAB binocdf

% Don' look at denormal numbers
kx = c > realmin;
err1 = err1(kx);
err2 = err2(kx);
err3 = err3(kx);
c = c(kx);
nk = nn(kx);

pltx = 1;
plt = [];
lg = {};
U = c;
figure(2)
subplot(2,1,1)
plt(pltx) = plot(nk,eps*(abs(log(U))),'k.');
lg{pltx} = 'Machine precision log CDF';
pltx = pltx + 1;
hold on

plt(pltx) = plot(nk,abs(err3./c),'.');
lg{pltx} = 'Relative error MATLAB binocdf vs Mathematica';
pltx = pltx + 1;


plt(pltx) = plot(nk,abs(err1./c),'.');
lg{pltx} = 'Relative error binomcdf-code vs Mathematica';
pltx = pltx + 1;

plt(pltx) = plot(nk,abs(err2./c),'.');
lg{pltx} = 'Relative error C-code vs Mathematica';
% pltx = pltx + 1;
set(gca,'YScale','log')
xlabel('k')
legend(plt, lg(:))

subplot(2,1,2)
plot(nk,(abs(Ch(1,kx)-c)./c)./(eps(abs(log(c)))),'.')
hold on
plot(nk,(abs(CM(1,kx)-c)./c)./(eps(abs(log(c)))),'.')
% plot(nk,(abs(Ch(2,kx)-c)./c)./(eps(abs(log(c)))),'.')
ylim([0,10])
legend({'BINOM-CDF','BINOCDF (MATLAB)'})
xlabel('k')
ylabel('relative error logCDF (ULP)')

% PMF evaluation
figure(3)
subplot(2,1,1)
U = Ph;
plot(nn,eps*(abs(log(U))),'k.')
hold on
load('pdata.mat');
P = P(:,1)';
plot(nn,abs((Ph - P)./P),'.')
plot(nn,abs((PM - P)./P),'.')
legend({'Log Prob','BINOM-PMF','BINOPDF (MATLAB)'})
xlabel('k')
set(gca,'YScale','log')

subplot(2,1,2)
plot(nn,(abs(Ph-P)./P)./(eps(abs(log(P)))),'.')
hold on
plot(nn,(abs(PM-P)./P)./(eps(abs(log(P)))),'.')
ylim([0,5])
legend({'BINOM-PMF','BINOPDF (MATLAB)'})
xlabel('k')
ylabel('relative error logPMF (ULP)')


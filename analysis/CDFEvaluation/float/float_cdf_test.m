%%
%  Evaluate accuracy of the Binomial CDF evaluation procedures in single
%  precision.
%  Compare against reference values in double precision from MATLAB's 
%  native binocdf (direct summation) and a MATLAB implementation of the 
%  binocdf.h C code (combination of summation and incomplete beta function 
%  evaluation).
%  
%  Note that the output of binom_cdf or binocdf_fast is more accurate than
%  the binocdf code in the tails of the distribution due to inaccuracies in
%  binopdf code, but this effect is not visible in single precision.
% 

addpath('../') % add double precision binom_cdf implementation
addpath('../../../mex')
if ~isfile('./../../mex/binoinv_fast.mexa64')
    mex ../../../mex/binoinv_fast.c -outdir ../../../mex -I../../../src/Serial
end

if ~isfile('./../../mex/binocdff_fast.mexa64')
    mex ../../../mex/binocdff_fast.c -outdir ../../../mex -I../../../src/Serial
end


%%
clear all
close all

N = 2^15;
% Make sure to choose p that does not suffer from floating point 
% representation error! Must take p = i/2^j for 0<=i<=2^j, i,j pos. int.
i = 120;
j = 10;
p = i/2^j
q = 1-p;

Ns = single(N);
ps = single(p);
qs = single(q);

nn = binoinv_fast(double(realmin('single')),N,p):binoinv_fast(1-double(eps(single(1))/2),N,p);

% CDF value array
Ch = zeros(3,length(nn));
CM = zeros(3,length(nn));
% PMF value array
Ph = zeros(2,length(nn));
PM = zeros(2,length(nn));

% Matlab implementation of the C-code CDF and PMF evaluation with single
% precision arguments
for ix = 1:length(nn)
    k = single(nn(ix));
    [a,b] = binom_cdff(k,Ns,ps,qs);   
    [aa,bb] = binom_cdf(k,Ns,ps,qs);
    if k > Ns*p
        Ch(1,ix) = b;
        Ch(2,ix) = bb;
    else
        Ch(1,ix) = a;
        Ch(2,ix) = aa;
    end
    Ph(1,ix) = binom_pmff(k,Ns,ps,qs);
    Ph(2,ix) = binom_pmf(k,Ns,ps,qs);
end
% Matlab native Binomial CDF evaluation (double precision)
CM(2,nn > N*p)  = binocdf(nn(nn >  N*p),N,p,'upper');
CM(2,nn <= N*p) = binocdf(nn(nn <= N*p),N,p);
% Matlab native Binomial CDF evaluation (single precision)
ns = single(nn);
CM(1,ns > Ns*ps)  = binocdf(ns(ns >  Ns*ps),Ns,ps,'upper');
CM(1,ns <= Ns*ps) = binocdf(ns(ns <= Ns*ps),Ns,ps);

% MEX version of C-code CDF evaluation
Ch(3,nn > N*p)  = binocdff_fast(single(nn(nn >  N*p)), Ns, ps, 1);
Ch(3,nn <= N*p) = binocdff_fast(single(nn(nn <= N*p)), Ns, ps);

% Matlab native PMF evaluation (double precision)
PM(2,:) = binopdf(nn,N,p);
% Matlab native PMF evaluation (single precision)
PM(1,:) = binopdf(ns,Ns,ps);

% Relative error PMF and CDF to MATLAB double precision values
errP = Ph - PM(2,:);
errP = errP./PM(2,:);

errC = Ch - CM(2,:);
errC = errC./CM(2,:);

%% Plotting
% CDF 
figure()
subplot(2,1,1)
plot(nn,abs(errC(1,:)),'.');
hold on
plot(nn,eps('single')*abs(log(single(CM(2,:)))),'k.')

legend({'relative error','\epsilon \times Log CDF'})
xlabel('k')
title('CDF evaluation')
set(gca,'YScale','log')

subplot(2,1,2)
plot(nn,abs(errC(1,:))./(eps(single(abs(log(CM(2,:)))))),'.')
xlabel('k')
ylabel('relative error logCDF (ULP)')


% PMF 
figure()
subplot(2,1,1)
plot(nn,abs(errP(1,:)),'.')
hold on
plot(nn,eps('single')*abs(log(single(PM(1,:)))),'k.')
legend({'relative error','\epsilon \times Log PMF'})
title('PMF evaluation')
xlabel('k')
set(gca,'YScale','log')

subplot(2,1,2)
plot(nn,abs(errP(1,:))./(eps(single(abs(log(PM(1,:)))))),'.')
xlabel('k')
ylabel('relative error logPMF (ULP)')

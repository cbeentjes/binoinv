%%
%
% Check what range the values ea = -lambda/a and eb = lambda/b take in the
% calculation of the incomplete beta function (betai). We use these values
% as the argument to rlog1 function in betaexp and betacf subroutines.
%
% It turns out -max(p,1-p)/2 <  ea <= 0 and
%                       0    <= eb < max(p, 1-p)
%

clear all
close all
N = 1e8;
p = 0.5;

% Note that we only call the beta function routine if 20 <= k <= N - 20 which
% constrains the values a & b can take.
b = 20:10:(N-20);
a = (N+1) - b;

k = b - 1;

y = p*ones(size(b));
x = 1 - y;

% In addition we only call the beta function routine whenever the following
% condition holds (N+1)p * 1/(2-p) + 1 < k + 1 < (N+1)p * 2/(1+p)
kx = and((b-1).*(2-y) > (N+1).*y, b.*(1 + y) < (N+1)*2.*y);

k = k(kx);

b = b(kx);
a = a(kx);
y = y(kx);
x = x(kx);


% Get correct lambda that is used in betai routines
lambda = zeros(size(b));
jx = a > b;
lambda(jx) = (N+1)*y(jx) - b(jx) + y(jx);
lambda(~jx) = (a(~jx) - (N+1)*x(~jx)) - x(~jx);
ix = lambda < 0;
lambda(ix) = - lambda(ix);
S = b(ix);
b(ix) = a(ix);
a(ix) = S;
S = x(ix);
x(ix) = y(ix);
y(ix) = S;

% Rlog1 arguments
ea = -lambda./a;
eb = lambda./b;

% Output extreme values
max(eb)
min(ea)

plot(k,eb)
hold on
plot(k,ea)
xlabel('k')
legend({'ea','eb'})

%% ERFCX argument
% Look at the argument z0 for erfcx in the betaexp routine (asymptotic
% expansion of the betainc function). 
%
% Note that z0 is always non-negative and given by z0 = sqrt(f) where
% f = a*rlog1(ea) + b*rlog(eb) under the condition that t=exp(-f) > 0. This
% automatically implies that f < 1074*log(2) < 750 so that we actually know
% that z0 < sqrt(750) < 30.
%
% We only call erfcx in betaexp, which we only call in the central region
% of the binomial distribution, given by the following conditions:
kx =  ( (a>100 & b >= a & 0.03*a >= lambda) | ...
         (b>100 & a >  b & 0.03*b >= lambda) );

% Apply restrictions
k = k(kx);
b = b(kx);
a = a(kx);
y = y(kx);
x = x(kx);
ea = ea(kx);
eb = eb(kx);

% Calculate and print max value for z0.
z0 = sqrt(a.*rlog1(ea) + b.*rlog1(eb));
max(z0)
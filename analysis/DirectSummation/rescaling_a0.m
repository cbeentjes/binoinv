%% 
% Let a0 = (1-p)^(-N) then we know that a0 is always a finite double/single
% precision variable when Np <= 10 or N(1-p) <= 10. 
% However, when Np > 10 and N(1-p) > 10 the value of a0 can cause overflow.
%
% Check for a given value of p for what value of N do we have either of the
% following conditions to be true:
% 1) P(K<10) > 2^-1074 (double) or 2^-149 (single)
% 2) P(K>N-10) > 2^-53 (double) or 2^-24 (single)
%
% If one of those conditions is satisfied then we are in the regime for
% which the algorithm would use direct summation. We check what the largest
% value of N is so that we can find a bound on the value of 
% a_0 = (1-p)^-N, which grows as N increases due to 0<p<1.
% 
% So we print for a given value of p the maximal N for which either
% condition holds and then the value of -N*log(1-p), which can be used to
% calculate a_0 via a_0 = exp(-N*log(1-p))
%

clear all
% p = 0.96015; % mu0 = 10 maximising value in double precision
% p = 0.86301; % mu0 = 10 maximising value in single precision
% p = 0.95131; % mu0 = 13 maximising value in double precision
% p = 0.8393279; % mu0 = 13 maximising value in single precision

p = 0.5;
%% Double precision
% Starting value to make sure that N*p > 10 and N*(1-p) > 10
mu0 = 10;
if p > 0.5
    N = ceil(mu0/(1-p));
else
    N = ceil(mu0/p);
end

M = [0,0];

% Note that for large Np and N(1-p) we only swap p and q if we are in the
% right tail of the distribution, i.e. if P(K>N-10) > 2^-53  
% So we must check whether it is possible for p and q to swap via flag2
flag1 = binocdf(9,N,p) > eps(0);
flag2 = binocdf(N-10,N,p,'upper') > eps(1)/2;

while flag1 || flag2  
    % Max value of a0 = (1-p)^-N depends on whether we can swap p & q for
    % the given combination of p and N
    if flag1 && flag2 % allowed both
        m = max(-N*log1p(-p), -N*log(p));
    elseif flag1
        m = -N*log1p(-p);
    else 
        m = -N*log(p);
    end
    
    if m > log(realmax) && M(1) <= log(realmax)
        N
        N*p
    end
    
    % Compare with current max
    if m > M(1)
        M(1) = m;
        M(2) = N;
    end        
    
    N = N + 1;
    
    flag1 = binocdf(9,N,p) > eps(0);
    flag2 = binocdf(N-10,N,p,'upper') > eps(1)/2;    
end

disp('Double precision')
disp(['N: ' num2str(M(2))])
disp(['log(a0): ' num2str(M(1))])

%% Single precision
mu0 = 10;
if p > 0.5
    N = ceil(mu0/(1-p));
else
    N = ceil(mu0/p);    
end

M = [0,0];

% Note that for large Np and N(1-p) we only swap p and q if we are in the
% right tail of the distribution, i.e. if P(K>N-10) > 2^-24  
% So we must check whether it is possible for p and q to swap via flag2
flag1 = binocdf(9,N,p) > eps(single(0));
flag2 = binocdf(N-10,N,p,'upper') > eps(single(1))/2;

while flag1 || flag2  
    % Max value of a0 = (1-p)^-N depends on whether we can swap p & q for
    % the given combination of p and N
    if flag1 && flag2 % allowed both
        m = max(-N*log1p(-p), -N*log(p));
    elseif flag1
        m = -N*log1p(-p);
    else 
        m = -N*log(p);
    end
    
    if m > log(realmax('single')) && M(1) <= log(realmax('single'))
        N
        N*p
    end    
    
    % Compare with current max
    if m > M(1)
        M(1) = m;
        M(2) = N;
    end    
    
    N = N + 1;
    
    flag1 = binocdf(9,N,p) > eps(single(0));
    flag2 = binocdf(N-10,N,p,'upper') > eps(single(1))/2;    
end

disp('Single precision')
disp(['N: ' num2str(M(2))])
disp(['log(a0): ' num2str(M(1))])
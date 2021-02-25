function p = binom_pmff(k,N,p,q)
% Matlab version of pmf evaluation routine used in binocdff.h header file
    if k == 0
        p = exp(N*log1p(-p));
        return
    end
    if k == N
        p = exp(N*log1p(-q));
        return
    end
    if k < 0 || k > N || k ~= round(k)
        p = single(0.0);
        return
    end
    
    lc = stirlerrf(N) - stirlerrf(k) - stirlerrf(N - k) ...
            - bd0f(k,N*p) - bd0f(N-k,N*q);
        
    if ( N > 2^63 )
        if k < single(0.6)*N
            lf = log(single(2*pi)) + log(k) + log1p(-k/N);
        else
            lf = log(single(2*pi)) + log(k) + log(N-k) - log(N);
        end
        p = exp(lc - single(0.5)*lf);
    
    else
        p = exp(lc)*sqrt(N/(single(2*pi)*k*(N-k)));        
    end
end

function s = stirlerrf(N)   
    
    sfe =  [
    single(0.0)         , single(0.081061467), ...
    single(0.041340696) , single(0.027677926), ...
    single(0.020790672) , single(0.016644691), ...
    single(0.013876129) , single(0.011896710),  ...
    single(0.010411265) , single(0.0092554622), ...
    single(0.0083305634), single(0.0075736755), ...
    single(0.0069428401), single(0.0064089942), ...
    single(0.0059513701), single(0.0055547336)];

    S0 = single(0.083333333);
    S1 = single(0.0027777778);

    if N < 16        
        s = sfe(N+1);
    else
        if N > 500
            s =  S0/N;
        else    
            s = (S0 -S1/(N*N))/N;        
        end
    end 

end

function s = bd0f(k, NP)
%   Coefficients from bd0Expansion.nb notebook.
    P0 = single(0.19999997);
    P1 = single(-0.11809043);
    Q1 = single(-1.3047534);
    Q2 = single(0.37667380);    
    
    if abs(k - NP) < single(0.5)*(k+NP)
        v = (k - NP)/(k + NP);
        s = (k - NP)*v;
        t = v*v;
        w = single(1/3) + t.*( P1.*t + P0) / ...
                  ((Q2.*t + Q1).*t + 1);       
    
        s = s + (single(2)*t*v*w)*k;
    else
        s = k*log(k/NP) + NP - k;
    end
end


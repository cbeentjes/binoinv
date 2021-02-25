function p = binom_pmf(k,N,p,q)
% Matlab version of pmf evaluation routine used in binocdf.h header file
    if k == 0
        p = exp(N*log1p(-p));
        return
    end
    if k == N
        p = exp(N*log(p));
        return
    end
    if k < 0 || k > N || k ~= round(k)
        p = 0.0;
        return
    end
    
    lc = stirlerr(N) - stirlerr(k) - stirlerr(N - k) ...
            - bd0(k,N*p) - bd0(N-k,N*q);
        
    if ( N > 2^63 )
        if k < 0.6*N
            lf = log(2*pi) + log(k) + log1p(-k/N);
        else
            lf = log(2*pi) + log(k) + log(N-k) - log(N);
        end
        p = exp(lc - 0.5*lf);
    
    else
        p = exp(lc)*sqrt(N/(2*pi*k*(N-k)));        
    end
end

function s = stirlerr(N)   
    
    sfe = [0.0                   , 0.081061466795327258219670264, ...
    0.041340695955409294093822081, 0.0276779256849983391487892927, ...
    0.020790672103765093111522771, 0.0166446911898211921631948653, ...
    0.013876128823070747998745727, 0.0118967099458917700950557241, ...
    0.010411265261972096497478567, 0.0092554621827127329177286366, ...
    0.008330563433362871256469318, 0.0075736754879518407949720242, ...
    0.006942840107209529865664152, 0.0064089941880042070684396310, ...
    0.005951370112758847735624416, 0.0055547335519628013710386899];

    S0 = 0.083333333333333333333;
    S1 = 0.00277777777777777777778;
    S2 = 0.00079365079365079365079365;
    S3 = 0.000595238095238095238095238;
    S4 = 0.0008417508417508417508417508;
    
    nn = N*N;
    if N < 16        
        s = sfe(N+1);
    else
        if N > 35
            if N > 80
                if N > 500
                    s =  (S0 - S1/nn)/N;
                else            
                    s = (S0 - (S1 - S2/nn)/nn)/N;        
                end
            else
                s = (S0 - (S1 - (S2 - S3/nn)/nn)/nn)/N;
            end
        else
            s = (S0 - (S1 - (S2 - (S3 - S4/nn)/nn)/nn)/nn)/N;        
        end
    end 

end

function s = bd0(k, NP)
%   Coefficients from bd0Expansion.nb notebook.
    P0 =  0.2000000000000000;
    P1 = -0.3669101082562791;
    P2 =  0.2056846557845579;
    P3 = -0.03394955490410171; 
    P4 =  0.00005007384892489343; 
    Q1 = -2.548836255566920;
    Q2 =  2.293465048754533;
    Q3 = -0.8464624069944372;
    Q4 =  0.1046656636182324; 
    
    if abs(k - NP) < 0.5*(k+NP)
        v = (k - NP)/(k + NP);
        s = (k - NP)*v;
        t = v*v;
        w = 1/3 + t.*((((P4*t + P3).*t + P2).*t + P1).*t + P0) / ...
                  ((((Q4*t + Q3).*t + Q2).*t + Q1).*t + 1);       
    
        s = s + (2*t*v*w)*k;
    else
        s = k*log(k/NP) + NP - k;
    end
end


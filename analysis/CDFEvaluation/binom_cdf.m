function [S0,S1] = binom_cdf(k,N,p,q)
% Matlab version of CDF evaluation routine used in binocdf.h header file
    if nargin < 4
        q = 1 - p;
    end
    if N ~= round(N) || N < 0 || isinf(N)
        S0 = NaN; S1 = NaN;
        return
    end
    if ~(p >= 0 && q >= 0)
        S0 = NaN; S1 = NaN;
        return
    end
    if isnan(k)
        S0 = NaN; S1 = NaN;
        return
    end
    if (k < 0)
        S0 = 0; S1 = 1;
        return
    end
    if (k >= N)
        S0 = 1; S1 = 0;
        return
    end
    if (p == 0)
        S0 = 1; S1 = 0;
        return
    end
    if (q == 0) 
        S0 = 0; S1 = 1;
        return
    end    
        
    if   k*(2-p) > (N+1)*p && ... 
        (k+1)*(1+p) < (N+1)*2*p && ...
         k > 20 && N - k > 20
        [S0,S1] = betai(q, p, N-k, k+1, N);
    else
        S = 0.0;

        if N < 2^53
            flag = k > N*p;
            if flag
                k = N - (k+1);
                n = ceil(k);
                a=p; p=q; q=a;
            else
                n = fix(k);
            end
            
            T = binom_pmf(n,N,p,q);            
            s = T;
            a = 1;
            t = T;
            for i = 0:49
                t = n*q*t;
                s = (N - n + 1)*p*s + t;
                a = (N - n + 1)*p*a; 
                n = n - 1;     
                if (a > 1e270)
                    a = a*1e-220;
                    S = S*1e-220;
                    T = T*1e-220;
                end
            end
            S = s/a;            
            
        else
            flag = (k > N*0.5);
            if flag
                k = N - k;
                a=p; p=q; q=a;
            end
            
            if k > N*p
                sgn = 1;
            else
                sgn = -1;
            end
            
            n = fix(k);
            
            flag = xor(flag, sgn > 0);
            
            n = n + flag*sgn;
            
            for i = 0:49
                T = binom_pmf(n + sgn*i, N, p, q);
                S = S + T;
            end
        end
        
        S0 = S;
        S1 = (0.5 - S) + 0.5;
        if flag
            a = S0; S0 = S1; S1 = a;
        end
    end      
end

function [S0,S1] = betai(x,y,a,b,apbm1)
    if a > b
        lambda = (apbm1*y - b) + y;
%         lambda = ((a+b-1)*y - b) + y;
%         lambda = (a+b)*y - b;
    else
        lambda = (a - apbm1*x) - x;
%         lambda = (a - (a+b-1)*x) - x;
%         lambda = a - (a+b)*x;
    end
    
    flag = lambda < 0;
    if flag
        S0 = b; b = a; a=S0;
        S0 = x; x = y; y=S0;
        lambda = abs(lambda);
    end   
    
    if ( (a>100 && b >= a && 0.03*a >= lambda) || ...
         (b>100 && a >  b && 0.03*b >= lambda) )
        S0 = betaexp(a,b,lambda);
    else
        S0 = betacf(x,y,a,b,lambda);
    end
    
    S1 = 0.5 + (0.5 - S0);
    if flag
        lambda = S0; S0 = S1; S1 = lambda;
    end
end

function S0 = betaexp(a,b,lambda)
    NUM = 20;
    E0 = 1.12837916709551;
    E1 = 0.353553390593274;
    
    A0 = zeros(NUM+1,1);
    B0 = zeros(NUM+1,1);
    C  = zeros(NUM+1,1);
    D  = zeros(NUM+1,1);
    
    eps = 100*1e-17;
    
    ea = -lambda/a;
    eb =  lambda/b;
    f = a*rlog1(ea) + b*rlog1(eb);

    t = exp(-f);
    if t == 0
        S0 = 0;
        return
    end
    
    z0 = sqrt(f);
    z = 0.5*(z0/E1);
    z2 = f + f;
    
    if (a < b)
        h = a/b;
        w0 = 1 / sqrt(a*(1+h));
        r1 = (b - a)/b;
    else
        h = b/a;
        w0 = 1 / sqrt(b*(1+h));
        r1 = (b - a)/a;
    end
    r0 = 1 / (1+h);
    
    A0(1) = (2/3)*r1;
    C(1) = -0.5*A0(1);
    D(1) = - C(1);
    % Custom erfcx implementation
%     j0 = (0.5/E0)*erfcx_implementation(z0);
    % Native MATLAB erfcx implementation
    j0 = (0.5/E0)*erfcx(z0);
    j1 = E1;
    sum = j0 + D(1)*w0*j1;
    
    s = 1;
    h2 = h*h;
    hn = 1;
    w = w0;
    znm1 = z;
    zn = z2;
    
    for n = 2:2:NUM
        hn = hn*h2;
        A0(n) = r0*2*(h*hn + 1)/(n + 2);
        s = s + hn;
        A0(n+1) = r1*2*s/(n+3);
        
        for i = n:n+1
            r = -0.5*(i+1);
            B0(1) = r*A0(1);
            for m = 2:i
                bsum = 0;
                for j=1:m-1
                    mmj = m - j;
                    bsum = bsum + (j*r - mmj)*A0(j)*B0(mmj);
                end           
                B0(m) = r*A0(m) + bsum/m;
            end
            C(i) = B0(i)/(i+1);      
        
            dsum = 0.0;
            for j = 1:i-1
                dsum = dsum + D(i-j)*C(j);
            end
            D(i) = -(dsum + C(i));
        end
        
        j0 = E1*znm1 + (n-1)*j0;
        j1 = E1*zn + n*j1;
        znm1 = znm1*z2;
        zn   = zn*z2;

        w = w*w0;
        t0 = D(n)*w*j0;
        w = w*w0;
        t1 = D(n+1)*w*j1;
        sum = sum + (t0 + t1);
        if abs(t0) + abs(t1) <= eps*sum
            break
        end
    end

    u = exp(-bcorr(a,b));
    S0 = E0*t*u*sum;    
end

function S0 = betacf(x,y,a,b,lambda)
    MAXIT = 1000;
    
    eps = 15*1e-17;
    bt = brcomp(a,b,lambda);

    if bt <= 0
        S0 = 0;
        return
    end
    
    c = 1 + lambda;
    c0 = b/a;
    c1 = 1 + 1/a;
    yp1 = y + 1;
    
    n = 0;
    p = 1;
    s = a + 1;    

    an0 = 0;
    bn0 = 1;
    an1 = 1;
    bn1 = c/c1;
    r1  = c1/c;
    
    while n <= MAXIT
        n = n + 1;
        t = n/a;
        w = n*(b-n)*x;
        e = a/s;
        alpha = (p*(p + c0)*e*e)*(w*x);
        e = (1+t)/(c1 + t + t);
        beta = n + w/s + e*(c + n*yp1);
        p = 1 + t;
        s = s + 2;
        
        t = alpha*an0 + beta*an1;
        an0 = an1; an1 = t;
        t = alpha*bn0 + beta*bn1;
        bn0 = bn1; bn1 = t;
        r0 = r1;
        r1 = an1/bn1;
        
        if abs(r1-r0) < eps*r1
            S0 = bt*r1;
            return
        end
        
        an0 = an0/bn1;
        bn0 = bn0/bn1;
        an1 = r1;
        bn1 = 1;
    end
    
    if n > MAXIT
        disp('Warning: reached max iterations in continued fractions expansion')
    end
    S0 = bt*r1;     
end

function B = brcomp(a,b,lambda)
    if a > b
        h = a/b;
        x0 = h/(1+h);
    else
        h = b/a;
        x0 = 1/(1+h);
    end
       
    ea = -lambda/a;
    eb =  lambda/b;
    u = rlog1(ea);
    v = rlog1(eb);
        
    z = exp(-(a*u + b*v)); 
    B = (1/sqrt(2*pi))*sqrt(b*x0)*z*exp(-bcorr(a,b));
end

function s = bcorr(a,b)
  
    C0 = 0.083333333333333333;
    C1 =-0.0027777777777760393;
    C2 = 0.00079365078287927436;
    C3 =-0.00059522072260950274;
    C4 = 0.00083168805235713158;


    if a > b
        h = b; b = a; a = h;
    end
    
    h = a/b;
    c = h/(1+h);
    x = 1/(1+h);
    x2 = x*x;
    
    S3 = 1 + (x + x2);
    S5 = 1 + (x + x2*S3);
    S7 = 1 + (x + x2*S5);
    S9 = 1 + (x + x2*S7);
    
    h = 1/b;
    t = h*h;
    
    w = (((C4*S9*t + C3*S7)*t + C2*S5)*t + C1*S3)*t + C0;
    
    w = w*(c/b);
    
    h = 1/a;
    t = h*h;
    
    s = ((((C4*t + C3)*t + C2)*t + C1)*t + C0)/a + w;    
    
end

function fp(x)
   fprintf("%.17e\n",x) 
end
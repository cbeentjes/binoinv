function [S0,S1] = binom_cdff(k,N,p,q)
% Matlab version of CDF evaluation routine used in binocdf.h header file
    if N ~= round(N) || N < 0 || isinf(N)
        S0 = NaN; S1 = NaN;
        return
    end
    if ~(p >= single(0.0) && q >= single(0.0))
        S0 = NaN; S1 = NaN;
        return
    end
    if isnan(k)
        S0 = NaN; S1 = NaN;
        return
    end
    if (k < single(0.0))
        S0 = single(0.0); S1 = single(1.0);
        return
    end
    if (k >= N)
        S0 = single(1.0); S1 = single(0.0);
        return
    end
    if (p == single(0.0))
        S0 = single(1.0); S1 = single(0.0);
        return
    end
    if (q == single(0.0)) 
        S0 = single(0.0); S1 = single(1.0);
        return
    end    
        
    if   k*(single(2)-p) > (N+single(1.0))*p && ... 
        (k+single(1.0))*(single(1.0)+p) < (N+1)*single(2)*p && ...
         k > single(20) && N - k > single(20)
        [S0,S1] = betai(q, p, N-k, k+single(1.0), N);
    else
        S = single(0.0);

        if N < single(2^24)
            flag = k > N*p;
            if flag
                k = N - (k+single(1.0));
                n = ceil(k);
                a=p; p=q; q=a;
            else
                n = fix(k);
            end
            
            T = binom_pmff(n,N,p,q);            
            s = T;
            a = 1;
            t = T;
            for i = 0:25
                t = n*q*t;
                s = (N - n + single(1.0))*p*s + t;
                a = (N - n + single(1.0))*p*a;
                n = n - single(1.0);   
                if (a > 1e20)
                    a = a*1e-10;
                    s = s*1e-10;
                    t = t*1e-10;
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
            
            for i = 0:19
                T = binom_pmff(n + sgn*i, N, p, q);
                S = S + T;
            end
        end
        
        S0 = S;
        S1 = (single(0.5) - S) + single(0.5);
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
    
    flag = lambda < single(0.0);
    if flag
        S0 = b; b = a; a=S0;
        S0 = x; x = y; y=S0;
        lambda = abs(lambda);
    end
    
    if ( (a>single(100) && b >= a && single(0.03)*a >= lambda) || ...
         (b>single(100) && a >  b && single(0.03)*b >= lambda) )
        S0 = betaexpf(a,b,lambda);
    else
        S0 = betacff(x,y,a,b,lambda);
    end
    
    S1 = single(0.5) + (single(0.5) - S0);
    if flag
        lambda = S0; S0 = S1; S1 = lambda;
    end
end

function S0 = betaexpf(a,b,lambda)
    NUM = 20;
    E0 = single(1.1283792);
    E1 = single(0.35355338);
    
    A0 = single(zeros(NUM+1,1));
    B0 = single(zeros(NUM+1,1));
    C  = single(zeros(NUM+1,1));
    D  = single(zeros(NUM+1,1));
    
    eps = single(100)*1e-7;
    
    ea = -lambda/a;
    eb =  lambda/b;
    f = a*rlog1f(ea) + b*rlog1f(eb);

    t = exp(-f);
    if t == single(0.0)
        S0 = single(0.0);
        return
    end
    
    z0 = sqrt(f);
    z = single(0.5)*(z0/E1);
    z2 = f + f;
    
    if (a < b)
        h = a/b;
        w0 = single(1.0) / sqrt(a*(single(1.0)+h));
        r1 = (b - a)/b;
    else
        h = b/a;
        w0 = single(1.0) / sqrt(b*(single(1.0)+h));
        r1 = (b - a)/a;
    end
    r0 = single(1.0) / (single(1.0)+h);
    
    A0(1) = (single(2)/single(3))*r1;
    C(1) = -single(0.5)*A0(1);
    D(1) = - C(1);
    j0 = (single(0.5)/E0)*erfcx(z0);
    j1 = E1;
    sum = j0 + D(1)*w0*j1;
    
    s = single(1.0);
    h2 = h*h;
    hn = single(1.0);
    w = w0;
    znm1 = z;
    zn = z2;
    
    for n = 2:2:NUM
        hn = hn*h2;
        A0(n) = r0*single(2)*(h*hn + single(1.0))/(n + single(2));
        s = s + hn;
        A0(n+1) = r1*single(2)*s/(n+single(3));
        
        for i = n:n+1
            r = -single(0.5)*(i+1);
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
        
            dsum = single(0.0);
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


    u = exp(-bcorrf(a,b));
    S0 = E0*t*u*sum;    
end

function S0 = betacff(x,y,a,b,lambda)
    MAXIT = 1000;
    
    eps = single(15.0)*1e-7;

    bt = brcompf(a,b,lambda);

    if bt <= single(0.0)
        S0 = single(0.0);
        return
    end
    
    c = single(1.0) + lambda;
    c0 = b/a;
    c1 = single(1.0) + single(1.0)/a;
    yp1 = y + single(1.0);
    
    n = single(0.0);
    p = single(1.0);
    s = a + single(1.0);    

    an0 = single(0.0);
    bn0 = single(1.0);
    an1 = single(1.0);
    bn1 = c/c1;
    r1  = c1/c;
    
    while n <= MAXIT
        n = n + single(1.0);
        t = n/a;
        w = n*(b-n)*x;
        e = a/s;
        alpha = (p*(p + c0)*e*e)*(w*x);
        e = (single(1.0)+t)/(c1 + t + t);
        beta = n + w/s + e*(c + n*yp1);
        p = single(1.0) + t;
        s = s + single(2.0);
        
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
        bn1 = single(1.0);
    end
    
    if n > MAXIT
        disp('Warning: reached max iterations in continued fractions expansion')
    end
    S0 = bt*r1;     
end

function B = brcompf(a,b,lambda)
    if a > b
        h = a/b;
        x0 = h/(single(1.0)+h);
    else
        h = b/a;
        x0 = single(1.0)/(single(1.0)+h);
    end
       
    ea = -lambda/a;
    eb =  lambda/b;
    u = rlog1f(ea);
    v = rlog1f(eb);
    
    z = exp(-(a*u + b*v)); 
    B = (1/sqrt(2*pi))*sqrt(b*x0)*z*exp(-bcorrf(a,b));
end

function s = bcorrf(a,b)

    C0 = single(0.083333333);
    C1 = single(-0.0027767253);

    if a > b
        h = b; b = a; a = h;
    end
    
    h = a/b;
    c = h/(single(1.0)+h);
    x = single(1.0)/(single(1.0)+h);
    
    S3 = single(1.0) + (x + x*x);

    
    h = single(1.0)/b;
    t = h*h;
    
    w = C1*S3*t + C0;
    w = w*(c/b);
    
    h = single(1.0)/a;
    t = h*h;
    
    s = (C1*t + C0)/a + w;    
end
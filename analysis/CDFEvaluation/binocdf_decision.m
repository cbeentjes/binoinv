%% 
% Type of different methods to evaluate binomial CDF. Value of y can be:
% -2    -   Direct summation (at start not p_{k-1}/p_k < 1/2)
% -1    -   Direct summation (p_{k-1}/p_k < 1/2)
%  1    -   Continued fraction incomplete beta function
%  2    -   Asymptotic expansion incomplete beta function

clear all
close all
% Binomial model parameters
N = 20000;
p = 0.4;
q = 1-p;

% Range of values for which the binomial pmf is (approx.) numerically > 0
N_l = floor(max(N*p - 40*sqrt(N*p*q),0));
N_h = ceil(min(N*p + 40*sqrt(N*p*q),N));

% Check what approximation method for the CDF is used for each value 
% within [N_l, N_h]
K = N_l:N_h;
Y = zeros(length(K),1);
ix = 1;
for k = K
    Y(ix) = bino_decision(k,N,p,q);    
    ix = ix + 1;
end

% Plot results
figure()
plot(K,Y)
hold on
% Add scaled version of the binomial PMF for reference
plot(K,2*binopdf(K,N,p)/binopdf(round(N*p),N,p))
xlabel('k')
xlim([N_l,N_h])
legend({'y','binomial pmf'})

%% Decision function for Binomial CDF

function y = bino_decision(k,N,p,q)
    % If not in tails of distribution evaluate incomplete beta function
    if (k*(2-p) > (N+1)*p && (k+1)*(1+p) < (N+1)*2*p) && (k > 20) && (N-k > 20)
        y = betai_decision(q,p,N-k,k+1,N);
    % If in tails of distribution use direct summation
    else
        % Swap
        if k > N*p
            k = N-(k+1);
            a=p; p=q; q=a;
        end
        % p_{k-1}/p_k >= 1/2
        if (k*(2-p)) >= (N+1)*p
            y = -2;
        % p_{k-1}/p_k < 1/2
        else
            y = -1;
        end
    end        
end


%% Decision function for incomplete beta function

function y = betai_decision(x,y,a,b,apbm1)
    if a > b
        lambda = (apbm1*y-b) + y;
    else
        lambda = (a - apbm1*x) - x;
    end
    if lambda < 0
        w=b; b=a; a=w;
        w=x; x=y; y=w;
        lambda = abs(lambda);
    end
    % Asymptotic expansion
    if (a > 100 && b >= a && lambda <= 0.03*a) || ...
       (b > 100 && a >  b && lambda <= 0.03*b)
        y = 2;
    % Continued fraction
    else
        y = 1;
    end
end
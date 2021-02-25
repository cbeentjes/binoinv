function y = rlog1f(x)
%
%   rlog1(x) = x - log(1+x)
%
%   Using minimax approximation for -0.618 < x < 1.618
%   Coefficients from rlog1Expansion.nb notebook.

    % Rational minimax coefficients 
    P0 = single(0.19999999);
    P1 = single(-0.11665086);
    Q1 = single(-1.2975472);
    Q2 = single(0.37142160);
    
    r = x./(x+single(2));
    t = r.*r;
    w = single(1/3) + t.*(P1.*t + P0) / ...
              ((Q2.*t + Q1).*t + single(1));
    y = t.*(x + single(2) - single(2)*r.*w);           
    
end
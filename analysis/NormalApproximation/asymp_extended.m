clear all

syms n p q k x t eps positive
syms w y z

% gam_l == log( sqrt(2*pi) * Gamma(x) )
% arg_l == log( integrand )

gam_4 = 0;

x = n+1;
gam_1 = (x-1/2)*log(x/n) - x + 1/(12*x) - 1/(360*x^3) + 1/(1260*x^5);   %  + O(1/n^7) error
gam_4 = gam_4 + (x-1/2)*log(n);

x = k;
gam_2 = (x-1/2)*log(x/n) - x + 1/(12*x) - 1/(360*x^3) + 1/(1260*x^5);   %  + O(1/n^7) error
gam_4 = gam_4 - (x-1/2)*log(n);

x = n-k+1;
gam_3 = (x-1/2)*log(x/n) - x + 1/(12*x) - 1/(360*x^3) + 1/(1260*x^5);   %  + O(1/n^7) error
gam_4 = gam_4 - (x-1/2)*log(n);

arg_l = gam_1 - gam_2 - gam_3 + gam_4 ...
      + (n-k)*log(t) + (k-1)*log(1-t) ...
      + (1/2)*(log(p)+log(1-p)-log(n));

% switch from t,k,n to z,y,eps, and expand in eps

arg_l = subs(arg_l,t,1-p+sqrt(p*(1-p)/n)*(z-y));
arg_l = subs(arg_l,k,n*p+sqrt(p*(1-p)*n)*y);
arg_l = subs(arg_l,n,eps^(-2));

arg_l = simplify(arg_l);

% Note that we can only take the Taylor expansion up to order 2*r - 1
% where r is the order of the Stirling series used in arg_l.

arg_l = simplify(taylor(arg_l,eps,'Order',8));

disp(' outcome of Taylor expansion')

% arg_l = arg_l / eps^2;

arg_l = simplify(arg_l)

term0 = simplify(subs(arg_l,eps,0))   % leading order term
arg = taylor(exp(arg_l-term0)-1,eps,'Order',8);  % error

% for each power of epsilon, decompose integrand error
% into powers of z, and then integrate by parts symbolically

for m=1:6
  arg  = diff(arg,'eps')/m;
  c    = simplify(subs(arg,eps,0));

  nmax = 0;
  while (diff(c,'z',nmax+1)~=0)
    nmax = nmax+1;
  end

  C{1} = subs(c,z,0);
  for n = 1:nmax-1
    C{n+1} = subs(diff(c,'z',n),z,0)/factorial(n);
  end
  C{nmax+1} = diff(c,'z',nmax)/factorial(nmax);

  f{m} = 0;
  for n = nmax:-1:1
    f{m} = f{m} - C{n+1}*y^(n-1);
    if (n>1)
      C{n-1} = C{n-1} + (n-1)*C{n+1};
    end
  end
end

f1 = simplify(f{1})
f2 = simplify(f{2})
f3 = simplify(f{3})
f4 = simplify(f{4})
f5 = simplify(f{5})
f6 = simplify(f{6})

% convert into asymptotic expansion for inverse

f1 = subs(f1,y,w);
f2 = subs(f2,y,w);
f3 = subs(f3,y,w);
f4 = subs(f4,y,w);
f5 = subs(f5,y,w);
f6 = subs(f6,y,w);

g0_1 = 1;
g0_2 = w;
g0_3 = 2*w*g0_2 + diff(g0_2,'w');
g0_4 = 3*w*g0_3 + diff(g0_3,'w');
g0_5 = 4*w*g0_4 + diff(g0_4,'w');
g0_6 = 5*w*g0_5 + diff(g0_5,'w');


g1   = - f1
g1_1 =            diff(g1,'w');
g1_2 =   w*g1_1 + diff(g1_1,'w');
g1_3 = 2*w*g1_2 + diff(g1_2,'w');
g1_4 = 3*w*g1_3 + diff(g1_3,'w');
g1_5 = 4*w*g1_4 + diff(g1_4,'w');


g2   = - g0_1*f2 - g0_2*f1^2/2 ...
       - g1_1*f1;
g2   = simplify(g2)
g2_1 =            diff(g2,'w');
g2_2 =   w*g2_1 + diff(g2_1,'w');
g2_3 = 2*w*g2_2 + diff(g2_2,'w');
g2_4 = 3*w*g2_3 + diff(g2_3,'w');


g3   = - g0_1*f3 - g0_2*f1*f2 - g0_3*f1^3/6 ...
       - g1_1*f2 - g1_2*f1^2/2 ...
       - g2_1*f1;
g3   = simplify(g3)
g3_1 =             diff(g3,'w');
g3_2 =   w*g3_1  + diff(g3_1,'w');
g3_3 = 2*w*g3_2  + diff(g3_2,'w');


g4   = - g0_1*f4 - g0_2*(f1*f3+f2^2/2) - g0_3*f1^2*f2/2 - g0_4*f1^4/24 ...
       - g1_1*f3 - g1_2*f1*f2 - g1_3*f1^3/6 ...
       - g2_1*f2 - g2_2*f1^2/2 ...
       - g3_1*f1;
g4   = simplify(g4)
g4_1 =            diff(g4,'w');
g4_2 =   w*g4_1 + diff(g4_1,'w');


g5   = - g0_1*f5 - g0_2*(f2*f3+f1*f4) - g0_3*(f1*f2^2/2 + f3*f1^2/2) ...
                                      - g0_4*f1^3*f2/6 - g0_5*f1^5/120 ...
       - g1_1*f4 - g1_2*(f1*f3+f2^2/2) - g1_3*f1^2*f2/2 - g1_4*f1^4/24 ...
       - g2_1*f3 - g2_2*f1*f2 - g2_3*f1^3/6  ...
       - g3_1*f2 - g3_2*f1^2/2 ...
       - g4_1*f1;
g5   = simplify(g5)
g5_1 =            diff(g5,'w');

g6   = - g0_1*f6 - g0_2*(f3^2/2 + f2*f4 + f1*f5) - g0_3*(f2^3/6 + f1*f2*f3 + f1^2*f4/2) ...
       - g0_4*(f1^2*f2^2/4 + f1^3*f3/6) - g0_5*(f1^4*f2/24) - g0_6*(f1^6/720) ...
       - g1_1*f5 - g1_2*(f2*f3+f1*f4) - g1_3*(f1*f2^2/2 + f3*f1^2/2) ...
                                      - g1_4*f1^3*f2/6 - g1_5*f1^5/120 ...
       - g2_1*f4 - g2_2*(f1*f3+f2^2/2) - g2_3*f1^2*f2/2 - g2_4*f1^4/24 ...
       - g3_1*f3 - g3_2*f1*f2 - g3_3*f1^3/6  ...
       - g4_1*f2 - g4_2*f1^2/2 ...
       - g5_1*f1;
g6   = simplify(g6)


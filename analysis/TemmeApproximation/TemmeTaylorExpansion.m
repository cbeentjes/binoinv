%%
%   Compare the accuracy of a Taylor expansion of the Temme uniform
%   asymptotic approximation for small eta. The evaluation of the Taylor
%   expansion is usually (much) cheaper than the use of Newton-Raphson
%   iteration to solve for the Temme approximation.
%   
%   The accuracy need not be exactly equal to that of the full
%   Newton-Raphson iteration because generally in the central region of the
%   distribution the delta-error bound is looser
%
%   Note eta = -norminv(u) / sqrt(N+1)
%
%   and for single precision we are interested in
%   CPU: expansion to order 3, |eta| < 6.0e-2*sqpq
%   GPU: expansion to order 9, |eta| < 5.8e-1*sqpq
%
%   for double precision we are interested in
%   CPU: expansion to order 5, |eta| < 1.25e-2*sqpq
%   GPU: expansion to order 9, |eta| < 5.48e-1*sqpq or 2.74*sqpq
%
clear all
close all

p = 1-1e-1;
% p = 0.0625;
% p = 0.125;
p=2^-10;
p = 1 - 2^-9;

% Order of Taylor expansion
n = 9;

% Use single precision or double precision
single_flag = true;

P = p;
if single_flag
    p = single(p);   
end

sqpq = sqrt(p - p*p);

% Set the limit of the Taylor expansion region, i.e. |eta| < y*sqpq
y = 5.8e-1;
% y = 6e-2;

% y = 5.48e-1;
% y = 2.74e-1;
% y = 7.5e-3;
% y = 1.25e-2;

y1 = max(-y,-sqrt(-2*log(p))/sqpq);
y2 = min(y,sqrt(-2*log1p(-p))/sqpq);

y1 = y1*(1-1e1*eps);
y2 = y2*(1-1e1*eps);

y1 = y1*sqpq;
y2 = y2*sqpq;

f = @(x) x.*log1p((P-x)./x) + (1-x).*log1p(-(P-x)./(1-x));
Y1 = fzero(@(x) f(x) + 0.5*y1^2,[P,1-eps]);
Y2 = fzero(@(x) f(x) + 0.5*y2^2,[1e-200,P]);

% exact answer
Nxi = 3e6;
xi_exact = linspace(Y1,Y2,Nxi);
F = f(xi_exact);
eta_exact = sqrt(-2*F);
eta_exact = sign(eta_exact).*sign(P-xi_exact).*eta_exact;

eta0 = eta_exact/sqpq;

%% Temme approximation via Taylor expansion
tic
pr = p-0.5;
if single_flag
    a(1) = single(1);
else
    a(1) = 1.0;
end
a(2) = pr*(1.0/3.0);                         
a(3) = (pr*pr-0.75)*(1.0/36.0);
a(4) = pr*(2.25-pr*pr)*(1.0/270.0);
a(5) = (-2.4375+pr*pr*(-13.5+pr*pr))*(1.0/4320.0);
a(6) = pr*(14.0625 + pr*pr*(19.5 + pr*pr))*(1.0/17010.0);
a(7) = -(13095. + pr*pr*(270756. + pr*pr*(143568. + 8896.*pr*pr)))*(1/348364800.);
a(8) = pr*(1287 + pr*pr*( 7692 + pr*pr*(1872 + 64*pr*pr)))/13063680;
a(9) = -(2049867 + pr*pr*(90331632 + pr*pr*(234206496 + pr*pr*(28628736 + ...
    146176.*pr*pr))))*(1./601974374400);
a(10) = pr*(5026455 + pr*pr*(67004496 + pr*pr*(90359712 + pr*pr*(5799168 - 71936*pr*pr))))*(1./387991296000);
a(11) = -(805522725 + pr*pr*(61005186612 + pr*pr*(369023296032 + ...
        pr*pr*(287140266624 + pr*pr*(9876089088 - 167812096*pr*pr)))))*(1/2224897287782400);
    
xi = 0;    
xi_eta = 0;
for i = 0:n-1
    xi = a(n-i) + eta0.*xi;
    xi_eta = (n-i)*a(n-i) + eta0.*xi_eta;
end
xi =  p - (p-p*p)*eta0.*xi;
xi_eta = xi_eta * -sqpq;
toc

%% Temme approximation via Newton iteration
if single_flag
    xi2 = single(p)*ones(size(xi));
else
    xi2 = p*ones(size(xi));
end
eta2 = zeros(size(xi));
tic
xi_eta2 = -sqpq;
K = 4;
for k = 0:K
    xi2 = xi2 - (eta2 - eta0*sqpq).*xi_eta2;
    xi2 = max(min(xi2,1 - 1/1e6),1/1e6);
    eta2 = sqrt(-2*(xi2.*log1p((p-xi2)./xi2) + (1-xi2).*log1p(-(p-xi2)./(1-xi2))));
    kx = p < xi2;
    eta2(kx) = - eta2(kx);
    xi_eta2 = -eta2 ./log1p((p-xi2)./((1-p)*xi2));
end
toc
% Note that the Newton iteration does not cope with eta0==0
xi2(eta0==0) = p;
xi_eta2(eta0==0) = -sqpq;

%% Plot
figure()
plot(eta0,abs(1-xi./xi_exact))
hold on
plot(eta0,abs(1-xi2./xi_exact),'.')
if single_flag
    yline(eps('single')/2,'r')
%     ylim([0,2e-7])    
else
    yline(eps,'r')
%     ylim([0,1e-10])
end
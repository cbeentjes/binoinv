clear all
close all

chebfun2eps eps

alpha = 0.15;
% alpha = 6e-2;

% p = chebfun2(@(p,s) p, [0,1,-alpha,alpha]);
% s = chebfun2(@(p,s) s, [0,1,-alpha,alpha]);

pr = chebfun2(@(pr,s) pr, [-0.5,0.5,-alpha,alpha]);
s = chebfun2(@(pr,s) s,  [-0.5,0.5,-alpha,alpha]);

a1 = 1;
a2 = pr*(1.0/3.0);
a3 = (pr.*pr-0.75)*(1.0/36.0);
a4 = pr.*(2.25-pr.*pr)*(1.0/270.0);
a5 = (-2.4375+pr.*pr.*(-13.5+pr.*pr))*(1.0/4320.0);
a6 = pr.*(14.0625 + pr.*pr.*(19.5 + pr.*pr))*(1.0/17010.0);
a7 = -(13095. + pr.*pr.*(270756. + pr.*pr.*(143568. + 8896.*pr.*pr)))*(1/348364800.);
a8 = pr.*(1287 + pr.*pr.*( 7692 + pr.*pr.*(1872 + 64*pr.*pr)))/13063680;
a9 = -(2049867 + pr.*pr.*(90331632 + pr.*pr.*(234206496 + pr.*pr.*(28628736 + ...
    146176.*pr.*pr))))*(1./601974374400);
a10 = pr.*(5026455 + pr.*pr.*(67004496 + pr.*pr.*(90359712 + pr.*pr.*(5799168 - 71936*pr.*pr))))*(1./387991296000);
a11 = -(805522725 + pr.*pr.*(61005186612 + pr.*pr.*(369023296032 + ...
        pr.*pr.*(287140266624 + pr.*pr.*(9876089088 - 167812096*pr.*pr)))))*(1/2224897287782400);
a12 = (657631629*pr + 15421266972*pr.^3 + 50785742880*pr.^5 + ...
 24381019008*pr.^7 + 453607680*pr.^9 - 5346304*pr.^11)/363159853056000;
a13 = (1/41650077227286528000).*(-1770433997445 - 204823343525400.*pr.^2 - ...
  2230890628651056*pr.^4 - 4442480886779136*pr.^6 - ...
  1378680477186816*pr.^8 - 13961811695616*pr.^10 + 21490970624*pr.^12);

xi = 0.5 + pr - (1/4 - pr.*pr).*(a1.*s + a2.*s.^2 + a3.*s.^3 + a4.*s.^4 + a5.*s.^5 + ...
    a6.*s.^6 + a7.*s.^7 + a8.*s.^8 + a9.*s.^9 + a10.*s.^10 + a11.*s.^11 + a12.*s.^12 + a13.*s.^13)

xi6 = 0.5 + pr - (1/4 - pr.*pr).*(a1.*s + a2.*s.^2 + a3.*s.^3 + a4.*s.^4 + a5.*s.^5 + ...
    a6.*s.^6);

A1 = @(pr) 1;
A2 = @(pr) pr*(1.0/3.0);
A3 = @(pr) (pr.*pr-0.75)*(1.0/36.0);
A4 = @(pr) pr.*(2.25-pr.*pr)*(1.0/270.0);
A5 = @(pr) (-2.4375+pr.*pr.*(-13.5+pr.*pr))*(1.0/4320.0);
A6 = @(pr) pr.*(14.0625 + pr.*pr.*(19.5 + pr.*pr))*(1.0/17010.0);
A7 = @(pr) -(13095. + pr.*pr.*(270756. + pr.*pr.*(143568. + 8896.*pr.*pr)))*(1/348364800.);
A8 = @(pr) pr.*(1287 + pr.*pr.*( 7692 + pr.*pr.*(1872 + 64*pr.*pr)))/13063680;
A9 = @(pr) -(2049867 + pr.*pr.*(90331632 + pr.*pr.*(234206496 + pr.*pr.*(28628736 + ...
    146176.*pr.*pr))))*(1./601974374400);
A10 = @(pr) pr.*(5026455 + pr.*pr.*(67004496 + pr.*pr.*(90359712 + pr.*pr.*(5799168 - 71936*pr.*pr))))*(1./387991296000);
A11 = @(pr) -(805522725 + pr.*pr.*(61005186612 + pr.*pr.*(369023296032 + ...
        pr.*pr.*(287140266624 + pr.*pr.*(9876089088 - 167812096*pr.*pr)))))*(1/2224897287782400);
A12 = @(pr) (657631629*pr + 15421266972*pr.^3 + 50785742880*pr.^5 + ...
 24381019008*pr.^7 + 453607680*pr.^9 - 5346304*pr.^11)/363159853056000;
A13 = @(pr) (1/41650077227286528000).*(-1770433997445 - 204823343525400.*pr.^2 - ...
  2230890628651056*pr.^4 - 4442480886779136*pr.^6 - ...
  1378680477186816*pr.^8 - 13961811695616*pr.^10 + 21490970624*pr.^12);


f_xi = @(pr,s) 0.5 + pr - (1/4 - pr.*pr).*(A1(pr).*s + ...
    A2(pr).*s.^2 + ...
    A3(pr).*s.^3 + ...
    A4(pr).*s.^4 + ...
    A5(pr).*s.^5 + ...
    A6(pr).*s.^6 + ...
    A7(pr).*s.^7 + ...
    A8(pr).*s.^8 + ...
    A9(pr).*s.^9 + ...
    A10(pr).*s.^10 + ...
    A11(pr).*s.^11 + ...
    A12(pr).*s.^12 + ...
    A13(pr).*s.^13);

res_xi = @(pr,s) (A3(pr).*s.^3 + ...
    A4(pr).*s.^4 + ...
    A5(pr).*s.^5 + ...
    A6(pr).*s.^6 + ...
    A7(pr).*s.^7 + ...
    A8(pr).*s.^8 + ...
    A9(pr).*s.^9 + ...
    A10(pr).*s.^10 + ...
    A11(pr).*s.^11 + ...
    A12(pr).*s.^12 + ...
    A13(pr).*s.^13);

XI = chebfun2(f_xi, [-0.5,0.5,-alpha,alpha],'eps',eps('single'))

R_XI = chebfun2(res_xi, [-0.5,0.5,-alpha,alpha],'eps',eps('single'))

[c,d,r] = cdr(xi);
[C,D,R] = cdr(XI);


plot(xi-XI)
% plot(xi - xi6);
xlabel('pr')
ylabel('s')
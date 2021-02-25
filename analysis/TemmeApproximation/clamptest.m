%%
%     References refer to 2020 paper by Gil, Segura & Temme
%     "ASYMPTOTIC INVERSION OF THE BINOMIAL AND NEGATIVE BINOMIAL
%     CUMULATIVE DISTRIBUTION FUNCTIONS"
%     DOI: 10.1553/etna_vol52s270
%
% The Temme approximation relies on the solution xi(eta) by solving the
% implicit equation (2.3) given by:
% -1/2*eta^2 = xi*log(p/xi) + (1-xi)*log((1-p)/(1-xi))
%
% Note that this can only be solved for xi when eta lies in the interval
% [ -sqrt(-2 log p), sqrt(-2 log(1-p)) ]
% 
% Here we check that if eta does NOT lie in this interval whether the
% clamped Newton iteration will return a value of xi either equal to 1/nu
% or 1-1/nu so that we get S = xi*nu is 1 or n-1, thus forcing the use of
% the direct summation algorithm.
%
% Observation 1)
% If the asymptotic approximation is used, i.e. |eta0| < sqrt(p(1-p))/100
% then we are sure that eta0 lies within the interval with a solution as we
% have the following |eta0| < sqrt(p(1-p)) < min(sqrt(-2 log(1-p)),sqrt(-2 log p))
%
clear all
clc

% Binomial parameters
p = 0.2;
N = ceil(10.01*max(1/p,1/(1-p)));

q    = 1-p;
sqpq = sqrt(p-p*p);
nu   = N+1;

% Try two extreme possibilities + boundaries of valid interval for eta0
U = [eps(0), ...
    normcdf(-sqrt(nu)*sqrt(-2*log1p(-p))), ...
    normcdf(sqrt(nu)*sqrt(-2*log(p))), ...
    1-eps(1)/2];
U(2) = U(2) - 50*eps(U(2));
U(3) = U(3) + 50*eps(U(3));
% Note that supplied U always satisfies 0 < U < 1
U = max(min(U,1-eps(1)/2),eps(0));

if N*p > 10 && N*q > 10
    for u = U
        % Starting value for eta given a U
        eta0 = -norminv(u) / sqrt(nu)

        % Valid interval for eta0
        disp('Valid interval for eta:')
        disp([-sqrt(-2*log(p)),sqrt(-2*log1p(-p))])

        % Check what happens if eta NOT in valid interval
        if eta0 < -sqrt(-2*log(p)) || eta0 > sqrt(-2*log1p(-p))
            disp('No solution to (2.3) possible, check what algorithm does:')
            if abs(eta0) < 1e-2*sqpq
                disp('asymptotic expansion')
                eta0 = eta0 / sqpq;
                pr = p -  0.5;
                % coeffs on page 273 below (2.7)
                a2 = pr*(1.0/3.0);
                a3 = (pr*pr - 0.75)*(1.0/36.0);
                a4 = pr*(2.25 - pr*pr)*(1.0/270.0);
                % expansion on page 273 in (2.7)
                xi = p - (p - p*p)*eta0*(1.0 + eta0*(a2 + eta0*(a3 + eta0*a4)));
                xi_eta = -sqpq*(1.0 + eta0*(2.0*a2 + eta0*(3.0*a3 + eta0*4.0*a4)));

                % coeffs and expansion on top of page 275
                a0 =  1.0/3.0    + xi*( -2.0/3.0);
                a1 =  1.0/36.0   + xi*(  5.0/36.0   + xi*(-5.0/36.0));
                a2 = -1.0/1620.0 + xi*(-21.0/1620.0 + xi*(69.0/1620.0 + xi*(-46.0/1620.0)));
                a3 =  7.0/6480.0 + xi*( -2.0/6480.0 + xi*(33.0/6480.0 + xi*(-62.0/6480.0 + xi*(31.0/6480.0))));
                eta0 = eta0*sqpq / sqrt(xi - xi*xi);
                eta1 = (-a0 + eta0*(a1 + eta0*(-a2 + eta0*a3))) / sqrt(xi - xi*xi);
            else
                disp('Newton method')
                flag = p > 0.5;
                if flag
                    a0=p; p=q; q=a0;
                    eta0 = -eta0;
                end
                
                xi = p;
                eta = 0;
                xi_eta = -sqpq;

                for k=0:4
                  xi = xi - (eta-eta0)*xi_eta;
                  xi = max(min(xi,1-1/nu),1/nu);   % clamp at two ends
                  eta = sqrt(-2*(xi*log1p((p-xi)/xi) + (1-xi)*log1p(-(p-xi)/(1-xi))));
                  eta = eta*sign(p-xi);
                  xi_eta = - eta / log1p((p-xi)/(q*xi));   % log1p used for accuracy
                end
                disp('Check LHS and RHS of (2.3)')
                disp([-0.5*eta0^2, xi*log(p/xi) + (1-xi)*log((1-p)/(1-xi))])

                eta1 = log(sqrt(xi-xi*xi)*eta0/(p-xi)) / eta0;
                
                if flag
                    xi = 1.0 - xi;
                    xi_eta = -xi_eta;
                    a0=p; p=q; q=a0;
                end
            end       
            x = xi*nu + eta1*xi_eta;
            % Check whether a valid solution gets returned. Only problematic
            % scenario is when 10 < x < N-10 is returned when (2.3) does not
            % have a solution.
            if 10 < x && x < N - 10
                disp([xi, eta1, xi_eta])
                disp('Seemingly valid solution found, but no solution exists!')
            else
                disp('Resort to direct summation. Succes.')
                disp(' ')
            end
        else
            disp('eta0 in valid interval, skip check.')
        end
    end
else
    disp('N*p or N*q is too small for Temme approximation to be used')
end

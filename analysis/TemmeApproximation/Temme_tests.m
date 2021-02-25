
clear all
close all
addpath('../CDFEvaluation')
addpath('../AS241_NormalInverseCDF')

addpath('../../mex')
if ~isfile('../../mexbinocdf_fast')
    mex ../../mex/binocdf_fast.c -outdir ../../mex -I../../src/Serial
end

addpath('../BRATIO')

% Nlist = [21 25 32 100 320 1000 3200];
% Nlist = [1e4, 1.5e4, 2e4, 2.5e4, 3e4, 6e4, 8e4];
% Nlist = [9e4, 1e5, 1.1e5, 1.2e5, 1.3e5];
% Nlist = [1.11e5, 1.12e5, 1.13e5];
% Nlist = [5e3,6e3,7e3,8e3,9e3,1e4];
Nlist = [9622];
% Nlist = [100 320 1000 3200 10000 50000];
% Nlist = [21 30 40];
% Nlist = [420 500 600 800];
% Nlist = [1e2,1e3,1e4,1e5,1e6,1e7];
% Nlist = [1e2,4e2,1e3,4e3,1e4,4e4,1e5,1e6,5e6];
% Nlist = [4e5,5e5,6e5,2e6];  
% Nlist = [1e6,2e6,4e6,8e6,1e7];
% Nlist = [2^25];
% Nlist = 165:180;
% Nlist = 1050:1080;
% Nlist = [150, 155, 160, 168, 169, 800, 1000, 1100, 2000, 2500];
% Nlist = 1900:100:2700;
% Nlist = 700:720;

CPU_flag = false;

single_flag = true;

num_p = 71; % Number of points for probability P to try
num_p = 151;

md = 0;
lg = {};

Err_data = cell(length(Nlist),num_p);
delta_data = cell(length(Nlist),num_p);
N_data = cell(length(Nlist),num_p);
p_data = cell(length(Nlist),num_p);

ix = 1;

tic
for N = Nlist

    ps   = [];
    errs = [];
    
    % Note that for CPU/GPU algorithm we switch from direct summation to
    % asymptotics when Np and N(1-p) > 10
    if 10/N >= 0.5
        ix = ix + 1;
        continue
    end
    plist = linspace(10/N, 1-10/N, num_p);
    
%     plist = linspace(10/N, 0.5);
    plist = linspace(0.9375,1-10/N,num_p);
%         plist = linspace(0.49,0.51,num_p);

    
    jx = 1;
    for p = plist

        sd = sqrt(N*p*(1-p));

        if single_flag
            N1 = max(10,floor(N*p + 1.4*norminv(eps(single(0)))*sd));
            N2 = min(floor(N*p + 1.4*norminv(1-eps(single(1))/2)*sd), N-10);
            N1 = double(N1);
            N2 = double(N2);
        else
            N1 = max(10,floor(N*p + 1.4*norminv(eps(0))*sd));
            N2 = min(floor(N*p + 1.4*norminv(1-eps(1)/2)*sd), N-10);
        end
            
        U  = [];
        X1 = [];
        X2 = [];
        X3 = [];
        R  = [];

        tmp = [];
        dN = 1;
        for n = N1:dN:N2
            u = binocdf_fast(n,N,p);
            v = 1 - u;

            % Technically we only need the Temme approximation for the CPU
            % algorithm when we are outside of the realm of the normal
            % asymptotic approximation. Note, however that for GPU
            % algorithm we always use the Temme approximation
            if CPU_flag && abs(norminv(u)) < 3
                continue
            end

            % Note that if u < realmin we will lose precision automatically
            % so ignore as this is the relevant effect we want to pick up
            % on here
            if single_flag
                if u < eps(single(0))
                  N1 = n + 1;
                  continue
                end
            else
                if u < realmin
                  N1 = n + 1;
                  continue
                end
            end

            % Similarly if u == 1 we can stop
            if u == 1
              N2 = n - 1;
              break
            end            
            
            % Rather than assuming that Q(u) = n we need to actually
            % calcuate r = Q(u) due to finite precision effects.            
            
            % To calculate the correct exact solution to Q(u) = r note that we
            % should use the complementary beta function for u close to 1,
            % otherwise we lose accuracy. We use
            % I_{1-p}(N-x+1,x)     +    I_{p}(x,N-x+1)            = 1 
            % betainc(1-p,N-x+1,x) + betainc(1-p,N-x+1,x,'upper') = 1
            % Note that the betainc function is not particularly accurate
            % for N large! Therefore we can try to hack around this issue
            % by wrapping the BRATIO function (TOMS708) from the R
            % standalone math library in a MEX function for large N
            if N < 1e4
                if u > 0.5
                  r = fzero(@(x) betainc(1-p,N-x+1,x,'upper') - v, [0,N]);
                else
                  r = fzero(@(x) betainc(1-p,N-x+1,x) - u, [n-2,n+2]);             
                end                   
            else
                if u > 0.5
                  r = fzero(@(x) bratio_wrapper(1-p,N-x+1,x,0) - v, [0,N]);              
                else
                  r = fzero(@(x) bratio_wrapper(1-p,N-x+1,x,1) - u, [n-2,n+2]);              
                end                   
            end
                       
            

            eta0 = -norminv_as241(u) / sqrt((N+1)*(p-p*p));
            % Zeroth order Temme approximation
            x1 = (N+1)*( p - (p-p*p)*eta0 );
            % First order Temme approximation            
            % FP32
            if single_flag 
                if CPU_flag
                    x2 = binominv(N,p,u,0,0);
                else
                    x2 = binominv(N,p,u,0,2);
                end
            % FP64
            else
                if CPU_flag
                    x2 = binominv(N,p,u,0,1);
                else
                    x2 = binominv(N,p,u,0,3);
                end
            end
            % Second order Temme approximation
            % FP32            
            if single_flag
                if CPU_flag
                    x3 = binominv(N,p,u,1,0);
                else
                    x3 = binominv(N,p,u,1,2);
                end                
            % FP64
            else
                if CPU_flag
                    x3 = binominv(N,p,u,1,1);
                else
                    x3 = binominv(N,p,u,1,3);
                end
            end

            U  = [U  u ];
            X1 = [X1 x1];
            X2 = [X2 x2];
            X3 = [X3 x3];

            R  = [R r];
            
            tmp = [tmp n];
            
        end
      
        % U-based delta forms (original Giles suggestion)
        
        if CPU_flag
            % Delta CPU
            % Note this delta needs to be valid for all u in [eps(0), 1-eps(1)/2]
            % when |abs(norminv(u)| > 3, i.e. the tails of the distribution.
            
            % Different forms to try for delta
%             delta = 6.1e-3 * (1-log(4*U.*(1-U)))./(N.*p.*(1-p));        
            delta = 2*(5.6e-3 * abs(1-2*p) + 4e-3/sqrt(N*p*(1-p))) * (1-0*log(4*U.*(1-U)))./(N.*p.*(1-p));        
        else
            % Delta GPU  
            % Note this delta needs to be valid for all u in [eps(0), 1-eps(1)/2]
            % whenever Q(u) > 10 and N - Q(u) > 10.    
            
            % Two different valid forms for delta. Note that this assumes
            % that Np > 15 and N(1-p) > 15.
%             delta = 0.0197*(1-log(4*U.*(1-U)))./(N.*p.*(1-p)); 
            if single_flag
                delta = 0.023* (abs(1-2*p) + 0.2/sqrt(N.*p.*(1-p))) *(1-4*log(4*U.*(1-U)))./(N.*p.*(1-p));                   
            else
                delta = 0.023* (abs(1-2*p) + 0.2/sqrt(N.*p.*(1-p))) *(1-log(4*U.*(1-U)))./(N.*p.*(1-p));   
            end
        end
        
        % X-based delta forms
        xi = X2./(N+1);
        delta = 0.0225*(1./((N+1)*(xi.*(1-xi))));

        % Can use a more p-specific form to make delta smaller, at cost of
        % a few more floating point operations (sqrt, abs)
%         if N > 500
%             delta = delta*(abs(1-2*p)*(1- 20/(sqrt(N))) + 20/(sqrt(N)));            
%         end
        
        % If delta is smaller than the relative machine precision of Q(u)
        % we actually don't really care about the precise value of delta
        % anymore. So cap delta from below at this value to machine
        % precision results.
        
        if single_flag
            delta = max(delta, eps(single(X3)));        
        else
            delta = max(delta, 5*eps(X3));
        end
        
        % If using FP64 GPU algorithm note that the Taylor expansion error
        % for larger eta0 must be corrected for
        if ~CPU_flag && ~single_flag
              eta0 = norminv(U)/sqrt((N+1)*p*(1-p));
              dx = and(abs(eta0) > 2.74e-1, abs(eta0) < 5.48e-1);              
              delta(dx) = max(4.55e-4,delta(dx));
        end              

        % Keep track of the max delta so that we can check whether it is
        % less than 1/2
        if ~isempty(delta)
            md = max(md,max(delta));      
        end

        if ~isempty(X3)
            ps   = [ps p];
            errs = [errs max(abs((X3 - R)./delta))];
        end
      
        % Check if max abs error is indeed smaller than delta
        [mrerr,mx] = max(abs((X3-R)./delta));
        if mrerr > 1
            disp('Error relative to delta exceeds 1!')
            fprintf(1,'N=%d, p=%f, relative err=%f \n',N,p,mrerr);
            delta(mx)
            abs(X3(mx)-R(mx))
            U(mx)
            0.58*sqrt((N+1)*(p-p*p))
            norminv(U(mx))/sqrt((N+1)*(p-p*p))
            % Check what the first value of eta0 is for which the relative
            % error is larger than 1.
            EE = abs((X3-R)./delta);
            EE(EE < 1) = inf;
            [~,EEx] = min(EE);
            norminv(U(EEx))/sqrt((N+1)*(p - p*p))                        
        end

        if p == 0.5
            fprintf(1,'N=%d, err=%g \n',N, max(abs((X3-R))));
        end

        Err_data{ix,jx} = X3 - R;
        delta_data{ix,jx} = delta;
        N_data{ix,jx} = tmp;
        p_data{ix,jx} = p;
        jx = jx + 1;
        
    end
    
    ix = ix + 1;

    figure(1)
    if ps
        plot(ps,errs); hold on
%         plot(N.*(1-ps).*ps,errs); hold on
        xlabel('p')
        ylabel('max error/bound')
        lg = [lg(:)' ['N=' num2str(N)]];
    end

end
toc
legend(lg)
ylim([0,1.1])

disp(['Max delta: ',num2str(md)])

%% Raw error plots

% close all
figure()
A = 0;
B = 0;

for jx = 1:length(Nlist)
% for jx = 6
    N = Nlist(jx);
    for ix = 1:num_p
        p = p_data{jx,ix};
        nplot = N_data{jx,ix};
        % Scale to xi = n/N to compare error for different values of N
        xi = nplot/N;

        plot(xi,abs(Err_data{jx,ix}),'.')
        
        hold on

%         plot(xi,delta_data{jx,ix},'.')
        % Plot the relative precision to check when the error is of the
        % same order
%         plot(xi,eps(nplot),'--')
%         plot(xi,-eps(nplot),'--')
        % Find the max error observed across the observations
        if max(abs(Err_data{jx,ix})) > A
            B = [N,p_data{jx,ix}];
            A = max(abs(Err_data{jx,ix}));
        end

    end
end
legend(lg)

%% Error relative to delta plots

% close all
figure()
% for jx = 1:length(Nlist)
for jx = length(Nlist)
    N = Nlist(jx);
    for ix = 1:num_p
        nplot = N_data{jx,ix};
        if isempty(nplot)
            continue
        end
                
        U = binocdf(N_data{jx,ix},N,p_data{jx,ix});
        E = abs(Err_data{jx,ix});
        p = p_data{jx,ix};
                
        ux = U<1;
        nplot = nplot(ux)/N;
        
        if CPU_flag
        % CPU
            % step 1
            uplot = (1-log(4*U(ux).*(1-U(ux))));

            % step 2
%             uplot = 0.0061*uplot;
            uplot = (0.0056*abs(1-2*p) + 0.004/sqrt(N*p*(1-p)))*uplot;

        else
        % GPU
            % step 1
            uplot = (1-4*log(4*U(ux).*(1-U(ux))));
            
            % step 2
%             uplot = 0.0197*uplot;        
%             uplot = 0.022*(abs(1-2*p) + 0.2/sqrt(N*p*(1-p)))*uplot;
            uplot = 0.023*(abs(1-2*p) + 0.2/sqrt(N*p*(1-p)))*uplot;

        end
        
        % step 3
        uplot = uplot/(N*p*(1-p));
        
        uplot = 2.25e-2*(1./(nplot*N) + 1./(N*(1-nplot))); 
%         uplot = uplot*(abs(1-2*p) + 20/(sqrt(N)));
%         uplot = uplot*(abs(1-2*p)*(1- 20/(sqrt(N))) + 20/(sqrt(N)));
        
        if ~CPU_flag && ~single_flag
              eta0 = norminv(U)/sqrt((N+1)*p*(1-p));
              dx = and(abs(eta0) > 2.74e-1, abs(eta0) < 5.48e-1);              
              uplot(dx) = max(4.55e-4,uplot(dx));
        end
        

        % Machine precision issues
        uplot = max(7*eps(nplot*N),uplot);
        
        % Check if the max delta is actually smaller than 1/2, which it
        % needs to be!
        if max(uplot) > 0.5
            N
            p
            max(uplot)
        end
        
        
        rerrplot = E(ux)./uplot;
        plot(nplot,rerrplot,'.')

        hold on
    end
end

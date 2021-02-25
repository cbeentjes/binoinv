%%
% MATLAB implementation of Temme approximation to binomial CDF inverse
%
% use fac to control whether to provide O(1) correction on top of the
% O(N) Temme uniform expansion.
%
% use domain_flag to control the extent of the Taylor approximation 
% to the Temme uniform expansion.
%
function x = binominv(N,p,U,fac,domain_flag)
  if nargin < 5
      domain_flag = 0;
  end

  q    = 1-p;
  sqpq = sqrt(p-p*p);
  nu   = N+1;
  
%   W = norminv_as241(U);
  W = norminv(U);
  eta0 = - W / sqrt(nu);  

  switch domain_flag
      case 3            % Double precision GPU-style
          Y = 3.6e-1;
      case 2            % Single precision GPU-style
          Y = 5.8e-1;
      case 1            % Double precision CPU-style
          Y = 1.25e-2;
      case 0            % Single precision CPU-style
          Y = 6.0e-2;
  end
  
  if domain_flag == 3
    Y = 5.48e-1;
    domain_flag = 2;
  end
              
  
  if (abs(eta0)<Y*sqpq)   % use asymptotics when eta near 0
    eta0 = eta0 / sqpq;
    pr = p-0.5;
    % first coefs on page 273 below (2.7)
    a2 = pr*(1.0/3.0);                         
    a3 = (pr*pr-0.75)*(1.0/36.0);
    a4 = pr*(2.25-pr*pr)*(1.0/270.0);
    a5 = (-2.4375+pr*pr*(-13.5+pr*pr))*(1.0/4320.0);
    % further coefs from TemmeCoefficients.nb file
    if domain_flag >= 2
        a6 = pr*(14.0625 + pr*pr*(19.5 + pr*pr))*(1.0/17010.0);
        a7 = -(13095. + pr*pr*(270756. + pr*pr*(143568. + 8896.*pr*pr)))*(1/348364800.);
        a8 = pr*(1287 + pr*pr*( 7692 + pr*pr*(1872 + 64*pr*pr)))/13063680;
        a9 = -(2049867 + pr*pr*(90331632 + pr*pr*(234206496 + pr*pr*(28628736 + ...
                146176.*pr*pr))))*(1./601974374400);
        a10 = pr*(5026455 + pr*pr*(67004496 + pr*pr*(90359712 + pr*pr*(5799168 - 71936*pr*pr))))*(1./387991296000);
        a11 = -(805522725 + pr*pr*(61005186612 + pr*pr*(369023296032 + ...
                    pr*pr*(287140266624 + pr*pr*(9876089088 - 167812096*pr*pr)))))*(1/2224897287782400);
                
        a12 = pr*(657631629 + pr*pr*(15421266972 + pr*pr*(50785742880 + ...
                    pr*pr*(24381019008 + pr*pr*(453607680 - pr*pr*5346304)))))*(1/363159853056000);
        a13 = -(1/41650077227286528000)*(1770433997445 + pr*pr*(204823343525400 + ...
                pr*pr*(2230890628651056 + pr*pr*(4442480886779136 + ...
                pr*pr*(1378680477186816 + pr*pr*(13961811695616 - pr*pr*21490970624))))));
    end

    % expansion on page 273 in (2.7)
     
    xi = 0;
    
    if domain_flag >= 4
        xi = a13 + eta0*xi;
        xi = a12 + eta0*xi;        
    end
    
    % FP64 GPU
    if domain_flag >= 3
        xi = a11 + eta0*xi;
        xi = a10 + eta0*xi;
    end    
    % FP32 GPU
    if domain_flag >= 2            
        xi = a9 + eta0*xi;
        xi = a8 + eta0*xi;
        xi = a7 + eta0*xi;
        xi = a6 + eta0*xi;           
    end
    % FP64 CPU
    if domain_flag >= 1
        xi = a5  + eta0*xi;
        xi = a4  + eta0*xi;
    end
    % FP32 CPU
    xi = a3  + eta0*xi;
    xi = a2  + eta0*xi;
    xi = 1.0 + eta0*xi;   
    xi_min_p =  - (p-p*p)*eta0*xi;
    
    xi = p + xi_min_p;   

    xi_eta = 1;
    eta1 = 0;    
    
    a0 = -2*pr/3;
    a1 = (3.75 - 11*pr*pr)*(1/36);
    a2 = pr*(-15.75 + 7*pr*pr)*(1/810);
    a3 =     (30.375 + pr*pr*(-2.25 + pr*pr*71))*(1/9720);
    if domain_flag >= 2
        a4 = pr*(-58.375 + pr*pr*(43 + pr*pr*38))*(1/15120);
        a5 = (104574.375 + pr*pr*(1104205.5 + pr*pr*(-822006 + pr*pr*147688)))*(1/391910400);
        a6 = pr*(-378783 + pr*pr*(-818028   + pr*pr*( 697392 + pr*pr*19904)))*(1/587865600);
        a7 = (2243295 + pr*pr*(65623068 + pr*pr*(25084368 + pr*pr*(-44529600 + ...
                pr*pr*(-2855936)))))*(1/75246796800);
        a8 = pr*(-3009884265 + pr*pr*(-23683185648 + pr*pr*(5930384544 + ...
                pr*pr*(7472169216 + pr*pr*322279168))))*(1/27935373312000);
        a9 = (1134850173975 + pr*pr*(63723825690372 + pr*pr*(190882709890272 + ... 
                pr*pr*(-114364903591296 + pr*pr*(-33113133925632 + pr*pr*(-469792173056))))))*(1/300361133850624000);        
    end
    
    if domain_flag >= 4
        eta1 = a9 + eta0*eta1;
        eta1 = a8 + eta0*eta1;  
    end
    
    % FP64 GPU
    if domain_flag >= 3      
        eta1 = a7 + eta0*eta1;
        eta1 = a6 + eta0*eta1;
    end        
    % FP32 GPU
    if domain_flag >= 2
        eta1 = a5 + eta0*eta1;                
        eta1 = a4 + eta0*eta1;
    end
    % FP64 CPU
    if domain_flag >= 1
        eta1 = a3 + eta0*eta1;
        eta1 = a2 + eta0*eta1;        
    end
    % FP32 CPU    
    eta1 = a1 + eta0*eta1;
    eta1 = a0 + eta0*eta1;      
  
  else                  % use Newton iteration when not, see (2.12)
    flag = (p > 0.5);
    if flag
        a=p; p=q; q=a;
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
    
    eta1 = log(sqrt(xi-xi*xi)*eta0/(p-xi)) / eta0;   % see (3.6)
                                                     
    if flag
        xi = 1 - xi;
        xi_eta = -xi_eta;
    end
    
    xi_min_p = xi - p;
  end    

  x = N*p + (xi_min_p*nu + (fac*eta1*xi_eta + p));
end

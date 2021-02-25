%% figure illustrating CDF inversion
close all
clear all

figure(1)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.2 0.8]; set(gcf,'pos',pos);

n = 20;
p = 0.25;

k = linspace(0,n);

plot(k,betainc(1-p,n-k+1,k),'k'); hold on
axis([0 10 0 1]);

k = 0:n;
y = binocdf(k,n,p);

ks = reshape([k; k],2*length(k),1);
ys = reshape([y; y],2*length(k),1);
plot(ks(2:end),ys(1:end-1),'k'); hold on
xlabel('$k$','Interpreter','latex'); 
ylabel('$u$','Interpreter','latex')
ax = gca; 
yticks(0:0.1:1)
set(ax,'TickLabelInterpreter','latex')

k  = 5.8;
u  = betainc(1-p,n-k+1,k);
x2 = 5.0;

[ax, ay] = dsxy2figxy(gca, [0 x2],[u u]);
annotation('arrow',ax,ay,'HeadWidth',6)
[ax, ay] = dsxy2figxy(gca, [x2 x2],[u 0]);
annotation('arrow',ax,ay,'HeadWidth',6)
[ax, ay] = dsxy2figxy(gca, [0 k],[u u]);
annotation('arrow',ax,ay,'HeadWidth',6)
[ax, ay] = dsxy2figxy(gca, [k k],[u 0]);
annotation('arrow',ax,ay,'HeadWidth',6)
[ax, ay] = dsxy2figxy(gca, [k x2],[0 0]);
annotation('arrow',ax,ay,'HeadWidth',6)

% print('-deps2','fig2.eps')

%
% MATLAB function obtained from
% http://www.mathworks.co.uk/matlabcentral/fileexchange/
%  30347-sigplot/content/sigplot/sigplot/BasicFunctions/dsxy2figxy.m
%

%% Accuracy Normal approximation for fixed n
%
% WARNING - Matlab betainc(1-p,n-ks+1,ks) is not accurate enough for n 
%           large (~1e5). To check accuracy of Normal approximation in that
%           case check NormalApproximationError.nb, Mathematica notebook
%           which is capable of running arbitrary precision.
%
% close all
clear all

figure(2)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 1.5]; set(gcf,'pos',pos);
figure(3)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.5]; set(gcf,'pos',pos);

N  = [40 1000];
p = 0.25;

W = linspace(-3,3,100);

for jx = 1:length(N)
    n = N(jx);
    ks = zeros(size(W));
    % Calculate k such that I_{1-p}(N-k+1,k) == U, where U = normcdf(W)
    for ix = 1:numel(ks)
        if W(ix) > 0
            ks(ix) = fzero(@(x) betainc(1-p,n-x+1,x,'upper')-normcdf(-W(ix)),[0,n]);            
        else
            ks(ix) = fzero(@(x) betainc(1-p,n-x+1,x)-normcdf(W(ix)),[0,n]);
        end
    end
    % Calculate the effective w from ks (rather than using W). This is 
    % slightly more stable for large values of n.
    w = norminv(betainc(1-p,n-ks+1,ks));
    % Compare to Normal approximations
    x1 = n*p + sqrt(n*p*(1-p)).*w;
    x2 = x1 + (0.5 + (1-2*p)*(w.^2-1)/6);
    x3 = x2 + ( (7*p*(1-p) - 1) - (p*(1-p) + 0.5)*w.*w).*w/(36*sqrt(n*p*(1-p)));
    x4 = x3 + ((2*p-1)*(p+1)*(p-2)*(3*w.^4+7*w.^2-16)/1620)./(n*p*(1-p));
    delta =   abs((w.^4/400 +    w.^2 /200  + 1/100 )*(1+p)*(2-p).* ...
        (abs(1-2*p) + 0.25./sqrt(n*p*(1-p))))./(n*p*(1-p));

    figure(2)
    subplot(3,2,jx); plot(W,ks-x2,'k'); 
    xlim([-3.5,3.5]); xticks(-3:1:3);
    xlabel('$w$','interpreter','latex'); ylabel('$Q-\widetilde{Q}_{N1}$','interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    title(['$n$ = ' num2str(n) ', $p$ = ' num2str(p)],'interpreter','latex')
    
    subplot(3,2,jx+2); plot(W,ks-x3,'k')
    xlim([-3.5,3.5]); xticks(-3:1:3);    
    xlabel('$w$','interpreter','latex'); ylabel('$Q-\widetilde{Q}_{N2}$','interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    y = ylim;
    ylim(1.5*y);
    
    subplot(3,2,jx+4); plot(W,ks-x4,'k')
    xlim([-3.5,3.5]); xticks(-3:1:3);    
    xlabel('$w$','interpreter','latex'); ylabel('$Q-\widetilde{Q}_{N3}$','interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    
    figure(3)
    subplot(1,2,jx);   plot(W,(ks-x3)./delta,'k')
    ylim([-1,1])
    xlim([-3.5,3.5]); xticks(-3:1:3);    
    xlabel('$w$','interpreter','latex'); ylabel('$(Q-\widetilde{Q}_{N2})/\delta$','interpreter','latex')
    title(['$n$ = ' num2str(n) ', $p$ = ' num2str(p)],'interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
end

% figure(2)
% print('-deps2','asymp1.eps');

% figure(3)
% print('-deps2','asymp3.eps');


%% Accuracy Normal approximation for fixed k
%
% WARNING - Matlab betainc(1-p,N-k+1,k) is not accurate enough for N 
%           large (~1e5). To check accuracy of Normal approximation in that
%           case check NormalApproximationError.nb, Mathematica notebook
%           which is capable of running arbitrary precision.
%
% close all
clear all

figure(4)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 1.5]; set(gcf,'pos',pos);
figure(5)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.5]; set(gcf,'pos',pos);

ks  = [10 100000];
p = 0.25;
W = linspace(-3,3,100);

for jx = 1:length(ks)
    k = ks(jx);
    N = zeros(size(W));
    % Calculate N such that I_{1-p}(N-k+1,k) == U, where U = normcdf(W)    
    for ix = 1:numel(N)
        if W(ix) > 0
            N(ix) = fzero(@(n) betainc(1-p,n-k+1,k,'upper')-normcdf(-W(ix)),[k-1,20*k/p]);            
        else
            N(ix) = fzero(@(n) betainc(1-p,n-k+1,k)-normcdf(W(ix)),[k-1,20*k/p]);
        end                
        % Note that N*p >= 10 and N*(1-p) >= 10 must hold if we are to use
        % the normal asymptotic approximation
        N(ix) = max([N(ix),10/p,10/(1-p)]);
    end
    % N must be integer!
    N = floor(N);
    % Calculate the effective w from N (rather than using W). This is 
    % slightly more stable for large values of N.
    w = norminv(betainc(1-p,N-k+1,k));
    
    % Compare to Normal approximations   
    x1 = N*p + sqrt(N*p*(1-p)).*w;
    x2 = x1 + (0.5 + (1-2*p)*(w.^2-1)/6);
    x3 = x2 + ( (7*p*(1-p) - 1) - (p*(1-p) + 0.5)*w.*w).*w./(36*sqrt(N*p*(1-p)));
    x4 = x3 + ((2*p-1)*(p+1)*(p-2)*(3*w.^4+7*w.^2-16)/1620)./(N*p*(1-p));
    delta =   abs((w.^4/400 +    w.^2 /200  + 1/100 )*(1+p)*(2-p).* ...
        (abs(1-2*p) + 0.25./sqrt(N*p*(1-p))))./(N*p*(1-p));    

    figure(4)
    subplot(3,2,jx); plot(N,k-x2,'k.')
    xlabel('$n$','interpreter','latex'); ylabel('$Q-\widetilde{Q}_{N1}$','interpreter','latex')
    title(['$k$ = ' num2str(k) ', $p$ = ' num2str(p)],'interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    subplot(3,2,jx+2); plot(N,k-x3,'k.')
    xlabel('$n$','interpreter','latex'); ylabel('$Q-\widetilde{Q}_{N2}$','interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    subplot(3,2,jx+4); plot(N,k-x4,'k.')
    xlabel('$n$','interpreter','latex'); ylabel('$Q-\widetilde{Q}_{N3}$','interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    
    figure(5)
    subplot(1,2,jx);   plot(N,(k-x3)./delta,'k.')
    ylim([-1,1])
    xlabel('$n$','interpreter','latex'); ylabel('$(Q-\widetilde{Q}_{N2})/\delta$','interpreter','latex')
    title(['$k$ = ' num2str(k) ', $p$ = ' num2str(p)],'interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
end

% figure(4)
% print('-deps2','asymp1.eps');
% 
% figure(5)
% print('-deps2','asymp3.eps');

%% Accuracy Normal approximation varying n
%
% WARNING - Matlab betainc(1-p,N-k+1,k) is not accurate enough for N 
%           large (~1e5). To check accuracy of Normal approximation in that
%           case check NormalApproximationError.nb, Mathematica notebook
%           which is capable of running arbitrary precision.
%
% close all
clear all

figure(6)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 1.0]; set(gcf,'pos',pos);

p = 0.25;
N = logspace(log10(max([10/p,10/(1-p)])),5,40);
W = linspace(-3,3,100);

for jx = 1:length(N)
    n = N(jx);
    % Calculate k such that I_{1-p}(N-k+1,k) == U, where U = normcdf(W)
    ks = zeros(size(W));
    for ix = 1:numel(ks)
        if W(ix) > 0
            ks(ix) = fzero(@(x) betainc(1-p,n-x+1,x,'upper')-normcdf(-W(ix)),[0,n]);            
        else
            ks(ix) = fzero(@(x) betainc(1-p,n-x+1,x)-normcdf(W(ix)),[0,n]);
        end        
    end
    % Calculate the effective w from ks (rather than using W). This is 
    % slightly more stable for large values of n.
    w = norminv(betainc(1-p,n-ks+1,ks));
    % Compare to Normal approximations
    x1 = n*p + sqrt(n*p*(1-p)).*w;
    x2 = x1 + (0.5 + (1-2*p)*(w.^2-1)/6);
    x3 = x2 + ( (7*p*(1-p) - 1) - (p*(1-p) + 0.5)*w.*w).*w/(36*sqrt(n*p*(1-p)));
    x4 = x3 + ((2*p-1)*(p+1)*(p-2)*(3*w.^4+7*w.^2-16)/1620)./(n*p*(1-p));
    delta =   abs((w.^4/400 +    w.^2 /200  + 1/100 )*(1+p)*(2-p).* ...
        (abs(1-2*p) + 0.25./sqrt(n*p*(1-p))))./(n*p*(1-p));
    
    % Calculate error relative to 'exact' solution ks (note ks is subject
    % to rounding error of betainc function!)
    e1(jx) = max(abs(x1-ks));
    e2(jx) = max(abs(x2-ks));
    e3(jx) = max(abs(x3-ks));
    e4(jx) = max(abs(x4-ks));
    re3(jx) = max(abs(x3-ks)./delta);
end


figure(6)
loglog(N,e2,'-.k',N,e3,'--k',N,e4,'xk')
% axis([10 1000 1e-4 1]) 
xlabel('$n$','interpreter','latex')
ylabel('maximum error','interpreter','latex')
legend('$||Q - \widetilde{Q}_{N1}||_\infty$','$||Q - \widetilde{Q}_{N2}||_\infty$','$||Q - \widetilde{Q}_{N3}||_\infty$','Location','Southwest',...
    'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

% print('-deps2','asymp2.eps');

figure(7)
semilogx(N,re3,'k')
ylim([0 1])
% axis([10 1000 0 1]) 
xlabel('$n$','interpreter','latex'); ylabel('maximum relative error','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

%% Accuracy Normal approximation varying k
%
% WARNING - Matlab betainc(1-p,N-k+1,k) is not accurate enough for N 
%           large (~1e5). To check accuracy of Normal approximation in that
%           case check NormalApproximationError.nb, Mathematica notebook
%           which is capable of running arbitrary precision.
%
% close all
clear all

figure(8)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 1.0]; set(gcf,'pos',pos);

% Can choose to load data straight from mathematica notebook
% (NormalApproximationError.nb). This data contains accurate combination of
% values k,M,p and w such that I_(1-p)(M-k+1,k) = normcdf(w).
load_mathematica = 1;

if load_mathematica
    load('mathematicadata.mat');
    ks = mathematicadata.xs;
    p = mathematicadata.p;
    W = mathematicadata.W';
    M = mathematicadata.M';
else
    p = 0.35;
    ks = logspace(-log10(p),-log10(p)+7,40);
    W = linspace(-3,3,10);
end

for jx = 1:length(ks)
    k = ks(jx);

    if load_mathematica
        N = M(:,jx);
        w = W(:,jx);
    else
        N = zeros(size(W));
        for ix = 1:numel(N)
            if W(ix) > 0
                N(ix) = fzero(@(n) betainc(1-p,n-k+1,k,'upper')-normcdf(-W(ix)),[k-1,20*k/p]);            
            else
                N(ix) = fzero(@(n) betainc(1-p,n-k+1,k)-normcdf(W(ix)),[k-1,20*k/p]);
            end                
            % Note that N*p >= 10 and N*(1-p) >= 10 must hold if we are to use
            % the normal asymptotic approximation
            N(ix) = max([N(ix),10/p,10/(1-p)]);          
        end
        % N must be integer!
        N = floor(N);
        % Calculate the effective w from N (rather than using W). This is 
        % slightly more stable for large values of N.
        w = norminv(betainc(1-p,N-k+1,k));
    end
    
    % Compare to Normal approximations   
    x1 = N*p + sqrt(N*p*(1-p)).*w;
    x2 = x1 + (0.5 + (1-2*p)*(w.^2-1)/6);
    x3 = x2 + ( (7*p*(1-p) - 1) - (p*(1-p) + 0.5)*w.*w).*w./(36*sqrt(N*p*(1-p)));
    x4 = x3 + ((2*p-1)*(p+1)*(p-2)*(3*w.^4+7*w.^2-16)/1620)./(N*p*(1-p));
    delta =   abs((w.^4/400 +    w.^2 /200  + 1/100 )*(1+p)*(2-p).* ...
        (abs(1-2*p) + 0.25./sqrt(N*p*(1-p))))./(N*p*(1-p));    
%     delta = max(delta,eps*(k));
    
    e1(jx) = max(abs(x1-k));
    e2(jx) = max(abs(x2-k));
    e3(jx) = max(abs(x3-k));
    e4(jx) = max(abs(x4-k));
    re3(jx) = max(abs(x3-k)./delta);
end

figure(8)
loglog(ks,e2,'-.k',ks,e3,'ok',ks,e4,'xk')
hold on
% Plot the infinite precision values from mathematica
if load_mathematica
    E = mathematicadata.Error;
    loglog(ks,E(2,:),'-.b',ks,E(3,:),'--b',ks,E(4,:),'xb')
    % Note that the accuracy of the approximations QN2 and QN3 seems to
    % deteriorate for large n, but this is due to the floating point error
    % in the calculation of ks, i.e. the answer is within floating point
    % precision of the correct answer.
    loglog(ks,eps(ks),'r') 
end 
hold off
% axis([10 1000 1e-4 1]) 
xlabel('$k$','interpreter','latex')
ylabel('maximum error','interpreter','latex')
legend('error 1','error 2','error 3','Location','Southwest','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

% print('-deps2','asymp2.eps');

figure(9)
semilogx(ks,re3,'k')
if load_mathematica
    hold on
    semilogx(ks,mathematicadata.RelativeError,'--')
    hold off
end
ylim([0 1])
% axis([10 1000 0 1]) 
xlabel('$k$','interpreter','latex'); ylabel('maximum relative error','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

%%

function varargout = dsxy2figxy(varargin)
% dsxy2figxy -- Transform point or position from data space 
% coordinates into normalized figure coordinates 
% Transforms [x y] or [x y width height] vectors from data space
% coordinates to normalized figure coordinates in order to locate
% annotation objects within a figure. These objects are: arrow, 
% doublearrow, textarrow, ellipse line, rectangle, textbox 
%
% Syntax:
%    [figx figy] = dsxy2figxy([x1 y1],[x2 y2])  % GCA is used
%    figpos      = dsxy2figxy([x1 y1 width height])
%    [figx figy] = dsxy2figxy(axes_handle, [x1 y1],[x2 y2])
%    figpos      = dsxy2figxy(axes_handle, [x1 y1 width height])
%
% Usage: Obtain a position on a plot in data space and  
%        apply this function to locate an annotation there, e.g., 
%   [axx axy] = ginput(2); (input is in data space)
%   [figx figy] = dsxy2figxy(gca, axx, axy);  (now in figure space)
%   har = annotation('textarrow',figx,figy); 
%   set(har,'String',['(' num2str(axx(2)) ',' num2str(axy(2)) ')']) 
%
%   Copyright 2006-2009 The MathWorks, Inc. 

% Obtain arguments (limited argument checking is done)
% Determine if axes handle is specified
if length(varargin{1}) == 1 && ishandle(varargin{1}) ...
                            && strcmp(get(varargin{1},'type'),'axes')	
	hAx = varargin{1};
	varargin = varargin(2:end); % Remove arg 1 (axes handle)
else
	hAx = gca;
end

% Remaining args are either two point locations or a position vector
if length(varargin) == 1        % Assume a 4-element position vector
	pos = varargin{1};
else
	[x,y] = deal(varargin{:});  % Assume two pairs (start, end points)
end

% Get limits
axun = get(hAx,'Units');
set(hAx,'Units','normalized');  % Make axes units normalized 
axpos = get(hAx,'Position');    % Get axes position
axlim = axis(hAx);              % Get the axis limits [xlim ylim (zlim)]
axwidth = diff(axlim(1:2));
axheight = diff(axlim(3:4));

% Transform from data space coordinates to normalized figure coordinates 
if exist('x','var')     % Transform a and return pair of points
	varargout{1} = (x - axlim(1)) * axpos(3) / axwidth + axpos(1);
	varargout{2} = (y - axlim(3)) * axpos(4) / axheight + axpos(2);
else                    % Transform and return a position rectangle
	pos(1) = (pos(1) - axlim(1)) / axwidth * axpos(3) + axpos(1);
	pos(2) = (pos(2) - axlim(3)) / axheight * axpos(4) + axpos(2);
	pos(3) = pos(3) * axpos(3) / axwidth;
	pos(4) = pos(4) * axpos(4 )/ axheight;
	varargout{1} = pos;
end

% Restore axes units
set(hAx,'Units',axun)

end
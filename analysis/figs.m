%
% first figure showing various regions in (x,p) space for fixed N
%

function figs()

close all
clear all

N = 1e4;

figure(1)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[0.8 0.8]; set(gcf,'pos',pos);

x1 = 10;
x2 = N-10;
x3 = 15/N;
x4 = 1-15/N;
xmin = 1;
xmax = N;
pmax = 1;
pmin = 1/N;


%
% bottom-up region
%

loglog(x1*[1 1],[x3 x4],'--k',...
       x2*[1 1],[x3 x4],'--k',...
       [x1 x2],x3*[1 1],'--k',...
       [x1 x2],x4*[1 1],'--k','LineWidth',1); 
hold on
axis square; 
axis([xmin xmax pmin pmax]);
xlabel('$k$','interpreter','latex'); ylabel('$p$','interpreter','latex'); 
text(1.5,0.05*pmax,'bottom-up','interpreter','latex')
text(5,1.5*pmin,'bottom-up','interpreter','latex')
text(10,2.5*pmin,'top-down','interpreter','latex')
h=text(2.5*x1,2.0*x3,'$|w|<3$','interpreter','latex');
set(h,'rotation',40);
h=text(7*x1,2.0*x3,'Temme','interpreter','latex');
h=text(1.5*x1,10.0*x3,'Temme','interpreter','latex');


set(gca,'TickLabelInterpreter','latex')

plot(binoinv(0.5,N,linspace(pmin,x3,1e3)),linspace(pmin,x3,1e3),'--k','LineWidth',1); 

%
% curves corresponding to different values of u
%

x = exp(linspace(log(1),log(xmax)-1e-10,100));


l_style = {'-.','-.',':',':','-','-'};
k = 1;
for u = [ normcdf(3) normcdf(-3) double(realmin('single')) -double(realmin('single')) 1e-300 -1e-300 ]

  if (u>0)   
    y = betaincinv(u,x+1,N-x,'upper');
  else
    y = betaincinv(-u,x+1,N-x,'lower');
  end
  if k > 2
      loglog(x,y,[l_style{k},'k'],'LineWidth',1);
  else
      loglog(x(x>=x1 & y>=x3),y(x>=x1 & y>=x3),[l_style{k},'k'],'LineWidth',1);
  end
  k = k + 1;
end
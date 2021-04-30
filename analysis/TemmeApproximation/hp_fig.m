%
% figure showing h_p(xi) function for various values of p
%
clear all
close all

xi = linspace(0,1,1e4);

figure()
pos = get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.6]; set(gcf,'pos',pos);

p = 1 - 0.125;
plot(xi,hp(xi,p),'k-.')

hold on

p = 0.5;
plot(xi,hp(xi,p),'k-')

p = 0.125;
plot(xi,hp(xi,p),'k--')

xlabel('$\xi$','interpreter','latex')
ylabel('$h_p(\xi)$','interpreter','latex')

ax = gca; 
set(ax,'TickLabelInterpreter','latex')
legend('$p=0.875$','$p=0.5$','$p=0.125$','interpreter','latex')

print('-deps2','approx1.eps')

%%
function h = hp(xi,p)
    h = sqrt(-2*(xi.*log1p((p-xi)./xi) + (1-xi).*log1p(-(p-xi)./(1-xi))));
    ix = p < xi;
    h(ix) = -h(ix); 
end
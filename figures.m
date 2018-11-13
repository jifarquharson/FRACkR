function [hFig] = figures(tt, P, time, KK, t, n, ke, w, kk, window_max)
%%   Plots three figures.
% (1): a multipanel figure showing the evolution of fracture porosity and permeability over time, the evolution of fracture permeability with porosity, and the evolution of equivalent permeability over time.
% (2): equivalent permeability evolution over time, contoured for different host rock permeabilities.
% (3): a schematic representation of the model geometry, indicating the randomly positioned fractures in terms of their depth.
hFig = figure(1);set(hFig, 'Position', [70 250 700 650]);
%%
% Plots porosity, permeability, equivalent permeability with time.
subplot (2,3,1);
semilogx((tt./86400),P); xlim([0 time]);
xlabel('Time [days]'), ylabel('Fracture porosity \phi_f');
subplot (2,3,2);
loglog((tt./86400), KK); xlim([0 time]);
xlabel('Time [days]'), ylabel('Fracture permeability k_f [m^2]');
title({['Permeability and porosity evolution after'];[num2str(round(t./86400), n),' days']});
subplot (2,3,3);
semilogy(P, KK);
xlabel('Fracture porosity \phi_f'), ylabel('Fracture permeability k_f [m^2]');
subplot (2,3,[4 6]);
semilogy((tt(:,1)./86400),ke); xlim([0 time]); ylim ([1*10^(-23) 1*10^(-10)]);
xlabel('Time [days]'), ylabel('Equivalent permeability k_e [m^2]');


hFig = figure(2);set(hFig, 'Position', [810 250 700 650]);
%%
% Plots equivalent permeability against time.
semilogy((tt(:,1)./86400), ke, 'k', 'Linewidth',1.1);
x = max((tt(:,1)./86400))*(4.1/5);
text(x, 2.5*10^-11, 'k_0 = 10^{-11}'); text(x, 2.5*10^-12, 'k_0 = 10^{-12}');
text(x, 2.5*10^-13, 'k_0 = 10^{-13}'); text(x, 2.5*10^-14, 'k_0 = 10^{-14}');
text(x, 2.5*10^-15, 'k_0 = 10^{-15}'); text(x, 2.5*10^-16, 'k_0 = 10^{-16}');
text(x, 2.5*10^-17, 'k_0 = 10^{-17}'); text(x, 2.5*10^-18, 'k_0 = 10^{-18}');
text(x, 2.5*10^-19, 'k_0 = 10^{-19}'); text(x, 2.5*10^-20, 'k_0 = 10^{-20}');
text(x, 2.5*10^-21, 'k_0 = 10^{-21}'); text(x, 2.5*10^-22, 'k_0 = 10^{-22}');
hold on, semilogy((tt(:,1)./86400), kk, 'b:', 'LineWidth',0.5);
xlabel('Time [days]'), ylabel('Equivalent permeability k_e [m^2]');
title(['Equivalent permeability evolution after ',num2str(t./86400),' days']);

cond = (tt(:,1)./86400);
fracs = repmat(w,[1,size(cond)]);
limit = 1.2*window_max; con_l = 0.25*max(cond); con_r = 0.5*max(cond);

hFig = figure(3); set(hFig, 'Position', [1550 250 300 650]);
%%
% Plots model geometry (i.e. fracture positions along conduit/dyke)
plot(cond,fracs);
set(gca,'YDir','Reverse');  ylim([0 limit]); set(gca,'xtick',[]);
ylabel('Depth [m]'); rectangle('Position',[con_l 0 con_r limit],'FaceColor',[1 0 0]);
text(0.5, 0.5, 'Conduit/dyke', 'rot',90, 'units','normalized');
text(0.125, 0.04, 'Host rock', 'rot',90, 'units','normalized');
text(0.875, 0.04, 'Host rock', 'rot',90, 'units','normalized');
title('Fracture positioning');

end


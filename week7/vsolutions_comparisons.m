% compare velocity solutions for different friction evolution laws
% Rishav Mallick, EOS, 2021

clear

% set up parameteric solutions
a = 1;
b = 100;
Vpl = 1;

% time vector
% t = linspace(0,10/a/Vpl,1e3)';
t = logspace(-4,log10(1/a),1e3)';

% velocity evo solutions
vlin = 1+exp(-a/.05.*t).*(b-1);
vpower2 = abs(tanh(atanh(sqrt(b)) + sqrt(Vpl)*a.*t).^2);
vlog = b./(b+exp(-Vpl*a.*t).*(1-b));

figure(1),clf
plot(t,vlin,'-','LineWidth',2,'Color',rgb('steelblue')), hold on
plot(t,vpower2,'-','LineWidth',2,'Color',rgb('lightgreen'))
plot(t,vlog,'-','LineWidth',2,'Color',rgb('purple'))
axis tight, grid on
legend('f = Av','f = Av^{1/2}','f = A log(v)')
xlabel('t'), ylabel('V(t)')
set(gca,'Fontsize',20,'YScale','lin','Xscale','lin','LineWidth',1)
% print(['vcomp_' num2str((Vpl)) 'vi_' num2str(b)],'-djpeg','-r300')

%% compute slip

slin = cumtrapz(t,vlin);
spower2 = cumtrapz(t,vpower2);
slog = cumtrapz(t,vlog);

figure(2),clf
plot(t,slin,'-','LineWidth',2,'Color',rgb('steelblue')), hold on
plot(t,spower2,'-','LineWidth',2,'Color',rgb('lightgreen'))
plot(t,slog,'-','LineWidth',2,'Color',rgb('purple'))
axis tight, grid on
legend('f = Av','f = Av^{1/2}','f = A log(v)')
xlabel('t'), ylabel('s(t)')
set(gca,'Fontsize',20,'YScale','lin','Xscale','lin','LineWidth',1)
% print(['scomp_' num2str((Vpl)) 'vi_' num2str(b)],'-djpeg','-r300')

%% compute stress

taulin = cumtrapz(t,a*(Vpl-vlin));taulin = taulin-min(taulin)+0.1;
taupower2 = cumtrapz(t,a*(Vpl-vpower2));taupower2 = taupower2 - min(taupower2)+0.1;
taulog = cumtrapz(t,a*(Vpl-vlog));taulog = taulog - min(taulog)+0.1;

figure(3),clf
plot(t,taulin,'-','LineWidth',2,'Color',rgb('steelblue')), hold on
plot(t,taupower2,'-','LineWidth',2,'Color',rgb('lightgreen'))
plot(t,taulog,'-','LineWidth',2,'Color',rgb('purple'))
axis tight, grid on
legend('f = Av','f = Av^{1/2}','f = A log(v)')
xlabel('t'), ylabel('\tau(t)')
set(gca,'Fontsize',20,'YScale','lin','Xscale','lin','LineWidth',1)

figure(4),clf
plot(vlin,taulin,'-','LineWidth',2,'Color',rgb('steelblue')), hold on
plot(vpower2,taupower2,'-','LineWidth',2,'Color',rgb('lightgreen'))
plot(vlog,taulog,'-','LineWidth',2,'Color',rgb('purple'))
axis tight, grid on
% plot power-law curves to compare
xplot = linspace(min(get(gca,'xlim')),max(get(gca,'xlim')),1e2)';
plot(xplot,.01*(xplot).^.1,'-.','LineWidth',2,'Color',rgb('purple'))
plot(xplot,.01*(xplot).^.5,'-.','LineWidth',2,'Color',rgb('lightgreen'))
plot(xplot,.01*(xplot).^1,'-.','LineWidth',2,'Color',rgb('steelblue'))
legend('f = Av','f = Av^{1/2}','f = A log(v)','\tau \propto v^{0.1}','\tau \propto v^{0.5}','\tau \propto v')

set(legend,'box','off','location','best')
xlabel('v'), ylabel('\tau')
ylim([1e-3 max(get(gca,'YLim'))])
set(gca,'Fontsize',20,'YScale','log','Xscale','log','LineWidth',1)

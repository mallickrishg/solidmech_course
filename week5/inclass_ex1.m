clear
addpath ~/Dropbox/scripts/utils/

% shear modulus
G = 30e3; %in MPa

% Stress kernels for fault slip
s12h=@(x2,x3,y2,y3,Wf) G*( ...
    -(x3-y3)./((x2-y2).^2+(x3-y3).^2)+(x3+y3)./((x2-y2).^2+(x3+y3).^2) ...
    +(x3-y3-Wf)./((x2-y2).^2+(x3-y3-Wf).^2)-(x3+y3+Wf)./((x2-y2).^2+(x3+y3+Wf).^2) ...
    )/2/pi;

s13h=@(x2,x3,y2,y3,Wf) G*( ...
    (x2-y2)./((x2-y2).^2+(x3-y3).^2)-(x2-y2)./((x2-y2).^2+(x3+y3).^2) ...
    -(x2-y2)./((x2-y2).^2+(x3-y3-Wf).^2)+(x2-y2)./((x2-y2).^2+(x3+y3+Wf).^2) ...
    )/2/pi;

% Displacement kernels for fault slip
u1h=@(x2,x3,y2,y3,W) ...
    (+atan((x3-y3)./(x2-y2))-atan((x3+y3)./(x2-y2)) ...
     -atan((x3-y3-W)./(x2-y2))+atan((x3+y3+W)./(x2-y2)) ...
    )/2/pi;

u1d=@(x2,x3,y2,y3,W,dip) ...
    (atan2((x3-y3),(x2-y2))-atan2((x3+y3),(x2-y2)) ...
    -atan2((x3-y3-W*sind(dip)),(x2-y2-W*cosd(dip)))+atan2((x3+y3+W*sind(dip)),(x2-y2-W*cosd(dip))) ...
    )/2/pi;
%% plot disp and strain for slip

W = 5e3;%m
% fault location
y2 = 0; 
y3 = 1e3;
dip = 50;
slip = 4;

x2 = linspace(-20e3,20e3,500);
x2long = x2';%linspace(-40e3,40e3,1000);
x3 = linspace(0e3,20e3,300);
[X2,X3] = meshgrid(x2,x3);

% compute displacement
u1 = u1h(X2,X3,y2,y3,W);
u1_dip = u1d(X2,X3,y2,y3,W,dip);

% compute stress
s12 = s12h(X2,X3,y2,y3,W);% in MPa
s13 = s13h(X2,X3,y2,y3,W);% in MPa

figure(1),clf
imagesc(x2./1e3,x3./1e3,u1), hold on
contour(x2./1e3,x3./1e3,u1,[-1:0.05:1].*0.5,'k-')
axis tight equal
xlabel('x_2 (km)'), ylabel('x_3 (km)')
cb=colorbar;
cb.Label.String = 'u_1 (m)';
set(gca,'YDir','reverse','FontSize',20)
colormap(bluewhitered(40))

figure(2),clf
plot(x2long./1e3,u1h(x2long,0,y2,y3,W),'-','LineWidth',2,'Color',rgb('orange')), hold on
plot(x2long./1e3,u1d(x2long,0,y2,y3,W,dip),'-','LineWidth',2,'Color',rgb('steelblue'))
axis tight, grid on
legend('vertical',['\delta = ' num2str(dip)])
xlabel('x_2 (km)'),ylabel('u_1 (m)')
set(gca,'FontSize',20,'LineWidth',1)

figure(10),clf
plot(x2long./1e3,slip.*u1h(x2long,0,y2,y3,W),'-','LineWidth',2,'Color',rgb('steelblue')), hold on
plot(x2long./1e3,0.*x2long,'k-','LineWidth',1)
axis tight, grid on
xlabel('x_2 (km)'),ylabel('u_1 (m)')
set(gca,'FontSize',15,'LineWidth',1)
print('ssfault_surfdis','-djpeg','-r200')
u1_surf = slip.*u1h(x2long,0,y2,y3,W);
writetable(table(x2long,u1_surf),'eq1_obs.dat')

figure(11),clf
plot(x2long./1e3,slip.*u1d(x2long,0,y2,y3,W,dip),'-','LineWidth',2,'Color',rgb('steelblue')), hold on
plot(x2long./1e3,0.*x2long,'k-','LineWidth',1)
axis tight, grid on
xlabel('x_2 (km)'),ylabel('u_1 (m)')
set(gca,'FontSize',15,'LineWidth',1)
print('ssdipfault_surfdis','-djpeg','-r200')
u1_surf = slip.*u1d(x2long,0,y2,y3,W,dip);
writetable(table(x2long,u1_surf),'eq2_obs.dat')

%% stresses
figure(3),clf
subplot(121)
imagesc(x2./1e3,x3./1e3,s12), hold on
contour(x2./1e3,x3./1e3,s12,[-5:.2:.5],'k-')
axis tight equal
xlabel('x_2 (km)'), ylabel('x_3 (km)')
cb=colorbar;
caxis([-1 1]*5)
cb.Label.String = '\sigma_{12} (MPa)';
set(gca,'YDir','reverse','FontSize',20)

subplot(122)
imagesc(x2./1e3,x3./1e3,s13), hold on
contour(x2./1e3,x3./1e3,s13,[-5:.2:.5],'k-')
axis tight equal
xlabel('x_2 (km)'), ylabel('x_3 (km)')
cb=colorbar;
caxis([-1 1]*5)
cb.Label.String = '\sigma_{13} (MPa)';
set(gca,'YDir','reverse','FontSize',20)
colormap bluewhitered(50)

%% compute max shear stress events
rng(42)
nev = 1e3;
x2ev = -10e3 + 20e3.*rand(nev,1);
x3ev = 1+19e3.*rand(nev,1);

dipev = 60+60.*rand(nev,1);
nvecev = [sind(dipev),-cosd(dipev)];

s12ev = s12h(x2ev,x3ev,y2,y3,W).*slip;
s13ev = s13h(x2ev,x3ev,y2,y3,W).*slip;

tvecev = s12ev.*nvecev(:,1) + s13ev.*nvecev(:,2);
ind = tvecev<1e-3;
tvecev(ind) = nan;

figure(14),clf
% contour(x2./1e3,x3./1e3,s12,[-5:.2:.5],'k-'), hold on
plot(zeros(100,1),linspace(0,20,100)','k-','LineWidth',2), hold on
scatter(x2ev./1e3,x3ev./1e3,30,tvecev,'filled','MarkerEdgeColor','none','LineWidth',.1)
axis tight equal, box on, grid on
xlabel('x_2 (km)'), ylabel('x_3 (km)')
cb=colorbar;
caxis(10.^[-1 1])
xlim([min(x2ev) max(x2ev)]./1e3)
cb.Label.String = '\Delta\sigma_{s} (MPa)';
set(gca,'YDir','reverse','FontSize',15,'ColorScale','log')
colormap parula
% print('aftershock_catalogue','-djpeg','-r200')

T = table(x2ev(~ind),x3ev(~ind),dipev(~ind),tvecev(~ind));
T.Properties.VariableNames = {'x2';'x3';'dip';'stressdrop_MPa'};
writetable(T,'aftershock_catalogue.dat')


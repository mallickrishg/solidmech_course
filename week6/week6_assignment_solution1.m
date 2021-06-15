% solutions to week 6 assignment part 1
% Rishav Mallick, EOS, 2021

clear

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

% Displacement kernels for a vertical fault
u1h=@(x2,x3,y2,y3,W) ...
    (+atan((x3-y3)./(x2-y2))-atan((x3+y3)./(x2-y2)) ...
     -atan((x3-y3-W)./(x2-y2))+atan((x3+y3+W)./(x2-y2)) ...
    )/2/pi;

% Displacement kernels for a dipping fault
u1d=@(x2,x3,y2,y3,W,dip) ...
    (atan2((x3-y3),(x2-y2))-atan2((x3+y3),(x2-y2)) ...
    -atan2((x3-y3-W*sind(dip)),(x2-y2-W*cosd(dip)))+atan2((x3+y3+W*sind(dip)),(x2-y2-W*cosd(dip))) ...
    )/2/pi;
%% load data

dat_eq = readtable('week6assignment/eq1_obs.dat');
aftershocks = readtable('week6assignment/aftershock_catalogue.dat');

figure(1),clf
subplot(311)
plot(dat_eq.x2long./1e3,dat_eq.u1_surf,'k.')
axis tight, grid on
xlabel('x_2 (km)'), ylabel('u_1 (m)')

subplot(3,1,[2,3])
scatter(aftershocks.x2./1e3,aftershocks.x3./1e3,50,aftershocks.stressdrop_MPa,'filled','MarkerEdgeColor','k')
axis tight equal, grid on, box on
caxis([0 1])
colormap parula
xlabel('x_2 (km)'), ylabel('x_3 (km)')
set(gca,'YDir','reverse')
%% create fault model
ss = [];
ss.N = 40;
y3edge = linspace(0,20e3,ss.N+1)';
W = diff(y3edge);
ss.W = W;
ss.y2c = zeros(ss.N,1);
ss.y3c = (W(1)/2:W:sum(W))';

% compute displace and traction kernels
Gd = zeros(length(dat_eq.x2long),ss.N);
K12 = zeros(length(aftershocks.x2),ss.N);
K13 = K12;

for i = 1:ss.N
    Gd(:,i) = u1h(dat_eq.x2long,0.*dat_eq.x2long,ss.y2c(i),ss.y3c(i)-ss.W(i)/2,ss.W(i));
    K12(:,i) = s12h(aftershocks.x2,aftershocks.x3,ss.y2c(i),ss.y3c(i)-ss.W(i)/2,ss.W(i));
    K13(:,i) = s13h(aftershocks.x2,aftershocks.x3,ss.y2c(i),ss.y3c(i)-ss.W(i)/2,ss.W(i));
end

%% do inverse problem

% linear inverse problems

% simple inversion
m1 = Gd\dat_eq.u1_surf;
% m1 = lsqnonneg(Gd,dat_eq.u1_surf);

figure(2),clf
subplot(211)
plot(dat_eq.x2long./1e3,dat_eq.u1_surf,'k.'), hold on
plot(dat_eq.x2long./1e3,Gd*m1,'r-','LineWidth',1)
axis tight, grid on
xlabel('x_2 (km)'), ylabel('u_1 (m)')
set(gca,'FontSize',20,'LineWidth',1)

subplot(212)
plot(m1,ss.y3c./1e3,'r-','LineWidth',1), hold on
axis tight, grid on
ylabel('x_3 (km)'), xlabel('Slip (m)')
set(gca,'YDir','reverse','FontSize',20,'LineWidth',1)

%% use aftershock information in inverse problem
% (1) use aftershocks to truncate Gd
ind = ss.y3c<.1e3 | ss.y3c>8e3;% remove these region from the design matrix
Gdm = Gd(:,~ind);
m2 = zeros(ss.N,1);
% m2(~ind) = Gdm\dat_eq.u1_surf;
m2(~ind) = lsqnonneg(Gdm,dat_eq.u1_surf);

figure(2),
subplot(211)
plot(dat_eq.x2long./1e3,Gd*m2,'-','LineWidth',1,'Color',rgb('forestgreen'))
subplot(212)
plot(m2,ss.y3c./1e3,'-','LineWidth',1,'Color',rgb('forestgreen'))

%% (2) try to match aftershock stress drops using K
n2 = repmat(sind(aftershocks.dip),1,ss.N);
n3 = repmat(-cosd(aftershocks.dip),1,ss.N);
Ktau = K12.*n2 + K13.*n3;

Gdm = [Gd(:,~ind);Ktau(:,~ind)];
m3 = zeros(ss.N,1);
% m3(~ind) = Gdm\[dat_eq.u1_surf;aftershocks.stressdrop_MPa];
m3(~ind) = lsqnonneg(Gdm,[dat_eq.u1_surf;aftershocks.stressdrop_MPa]);
stressdrop_pred = Ktau*m3;

figure(2),
subplot(211)
plot(dat_eq.x2long./1e3,Gd*m3,'-','LineWidth',1,'Color',rgb('royalblue'))
subplot(212)
plot(m3,ss.y3c./1e3,'-','LineWidth',1,'Color',rgb('royalblue'))

% figure(3),clf
% scatter(aftershocks.x2./1e3,aftershocks.x3./1e3,50,Ktau*m3,'filled','MarkerEdgeColor','k')
% axis tight equal, grid on, box on
% caxis([0 1])
% colormap parula
% xlabel('x_2 (km)'), ylabel('x_3 (km)')
% set(gca,'YDir','reverse')

%% nonlinear inverse problem
fun = @(m,x) m(1).*u1h(x,0.*x,0,m(2),m(3));
minit = [5,3e3,2e3];
mfit = nlinfit(dat_eq.x2long,dat_eq.u1_surf,fun,minit);
m4 = zeros(ss.N,1);
m4(ss.y3c-ss.W/2>=mfit(2) & ss.y3c+ss.W/2<=mfit(2)+mfit(3)) = mfit(1);

figure(2)
subplot(211)
plot(dat_eq.x2long./1e3,fun(mfit,dat_eq.x2long),'--','LineWidth',1,'Color',rgb('gray'))
subplot(212)
plot(m4,ss.y3c./1e3,'-','LineWidth',3,'Color',rgb('gray'))

%% plot dipping fault comparison
dat_eq2 = readtable('week6assignment/eq2_obs.dat');

dipvec = [30:5:90];
cspec = jet(length(dipvec));
lgdvec = {};

figure(3),clf
plot(dat_eq2.x2long./1e3,dat_eq2.u1_surf,'ks'), hold on
lgdvec{1} = 'observations';
for i = 1:length(dipvec)
    plot(dat_eq.x2long./1e3,mfit(1).*u1d(dat_eq2.x2long,0.*dat_eq2.x2long,0,mfit(2),mfit(3),dipvec(i)),'-','LineWidth',1,'Color',cspec(i,:))
    lgdvec{i+1} = ['\delta = ' num2str(dipvec(i))];
end
axis tight, box on, grid on
xlabel('x_2 (km)'), ylabel('u_1 (m)')
legend(lgdvec,'location','best')
set(gca,'FontSize',20,'LineWidth',1)

%% displacement kernels plot
figure(10),clf
pcolor(dat_eq.x2long./1e3,ss.y3c./1e3,abs(Gd)'), shading flat
hold on
contour(dat_eq.x2long./1e3,ss.y3c./1e3,abs(Gd)',logspace(-3,0,30),'k')
colormap(parula(150))
caxis(10.^[-4 -1]) 
xlim([0,max(dat_eq.x2long./1e3)])
xlabel('Observation point x_2 (km)'), ylabel('Fault patch x_3 (km)')
set(gca,'ColorScale','log','YDir','reverse','FontSize',20,'LineWidth',2)


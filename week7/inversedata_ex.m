% Inversion exercise
% Rishav Mallick, EOS, 2021

clear
addpath ~/Dropbox/scripts/utils

eq_dat = readtable('assignment_week7/eq_displacement.dat');
gps = readtable('assignment_week7/gps_timeseries.dat');

% eq_dat = readtable('eq_displacement.dat');
% gps = readtable('gps_timeseries.dat');

% extract important variables
ox = eq_dat.ox;
ueq = eq_dat.u1_eq;
t = gps.tvals;
usta = gps{:,2:end};

%% define greens functions
% shear modulus
G = 30e3; %in MPa

% Stress kernels for fault slip
s12h=@(x2,x3,y2,y3,Wf) G*( ...
    -(x3-y3)./((x2-y2).^2+(x3-y3).^2)+(x3+y3)./((x2-y2).^2+(x3+y3).^2) ...
    +(x3-y3-Wf)./((x2-y2).^2+(x3-y3-Wf).^2)-(x3+y3+Wf)./((x2-y2).^2+(x3+y3+Wf).^2) ...
    )/2/pi;

% Displacement kernels for a vertical fault
u1h=@(x2,x3,y2,y3,W) ...
    (+atan((x3-y3)./(x2-y2))-atan((x3+y3)./(x2-y2)) ...
     -atan((x3-y3-W)./(x2-y2))+atan((x3+y3+W)./(x2-y2)) ...
    )/2/pi;
%% create fault model
ss = [];
ss.N = 100;
y3edge = linspace(0,20e3,ss.N+1)';
W = diff(y3edge);
ss.W = W;
ss.y2c = zeros(ss.N,1);
ss.y3c = (W(1)/2:W:sum(W))';

% compute displace and traction kernels
Gd = zeros(length(ox),ss.N);
K12 = zeros(ss.N,ss.N);

for i = 1:ss.N
    Gd(:,i) = u1h(ox,0.*ox,ss.y2c(i),ss.y3c(i)-ss.W(i)/2,ss.W(i));
    K12(:,i) = s12h(ss.y2c,ss.y3c,ss.y2c(i),ss.y3c(i)-ss.W(i)/2,ss.W(i));
end
%% invert for earthquake slip
lvec = 10.^[-4:.1:2]';
meqvec = zeros(ss.N,length(lvec));
resvec = zeros(length(lvec),1);

% define smoothing matrix
% Lsm = eye(ss.N);
Lsm = K12;
% Lsm = compute_laplacian1d(ss.N);

for i = 1:length(lvec)
%     mvec(:,i) = inv(Gd'*Gd + lvec(i)^2*(Lsm'*Lsm))*Gd'*ueq;
    meqvec(:,i) = lsqnonneg([Gd;lvec(i).*Lsm],[ueq;zeros(ss.N,1)]);
    resvec(i) = (ueq-Gd*meqvec(:,i))'*(ueq-Gd*meqvec(:,i));
end

% choose models that are acceptable
mst = 1; mend = 5;

figure(10),clf
subplot(311)
semilogx(lvec,resvec,'-v','LineWidth',2), hold on
axis tight, grid on
ylabel('\epsilon^t \epsilon'), xlabel('\lambda')
set(gca,'FontSize',20,'LineWidth',1)
  
% plot slip results
figure(11),clf
subplot(2,1,1)
plot(ox./1e3,ueq,'ro'), hold on
plot(ox./1e3,Gd*meqvec(:,mst:mend),'-','LineWidth',1)
axis tight, grid on
xlabel('x_2 (km)'), ylabel('u_1 (m)')
set(gca,'FontSize',20,'LineWidth',1)

subplot(2,1,2)
plot(meqvec(:,mst:mend),ss.y3c./1e3,'-','LineWidth',1)
axis tight, grid on
xlabel('Slip (m)'), ylabel('x_3 (km)')
set(gca,'FontSize',20,'LineWidth',1,'YDir','reverse')

% plot stress change on fault
figure(12),clf    
plot(K12*meqvec(:,mst:mend),ss.y3c./1e3,'-','LineWidth',1), hold on
plot([0 0],get(gca,'YLim'),'k-','LineWidth',2)
axis tight, grid on
xlabel('\Delta\tau (MPa)'), ylabel('x_3 (km)')
set(gca,'FontSize',20,'LineWidth',1,'YDir','reverse')

%% interseismic timeseries
mintervec = zeros(ss.N,length(lvec));
rintervec = nan(length(lvec),1);

% Fault slip rate
Vpl = 1e-9*86400;% m/day

% time window
tind = t<-30;
tvec = t(tind);
ustamod = usta(tind,:)-usta(find(tind,1),:);

% design matrix
Gvel = [ones(length(tvec),1) tvec];
% estimate interseismic velocity for each station
vinter = zeros(length(ox),1);
for i = 1:length(ox)
    mvel = Gvel\ustamod(:,i);
    vinter(i) = mvel(2);
end
% remove arc-tangent
vinter = vinter - Vpl/pi*atan2(ox,20e3);

% invert for interseismic sliprate
for i = 1:length(lvec)
    mintervec(:,i) = lsqnonneg([Gd;lvec(i).*Lsm],[vinter;zeros(ss.N,1)]);
    rintervec(i) = (vinter-Gd*mintervec(:,i))'*(vinter-Gd*mintervec(:,i));
end

figure(10)
subplot(312),cla
semilogx(lvec,rintervec,'-v','LineWidth',2), hold on
axis tight, grid on
ylabel('\epsilon^t \epsilon'), xlabel('\lambda')
set(gca,'FontSize',20,'LineWidth',1)

mst = 1; mend = 5;
% plot interseismic results
figure(8),clf
subplot(2,1,1)
plot(ox./1e3,vinter+Vpl/pi*atan2(ox,20e3),'ro'), hold on
plot(ox./1e3,Gd*mintervec(:,mst:mend)+Vpl/pi*atan2(ox,20e3),'-','LineWidth',1)
axis tight, grid on
xlabel('x_2 (km)'), ylabel('v_1 (m/day)')
set(gca,'FontSize',20,'LineWidth',1)

subplot(2,1,2)
plot(mintervec(:,mst:mend)./Vpl,ss.y3c./1e3,'-','LineWidth',1)
axis tight, grid on
xlabel('Normalized Sliprate'), ylabel('x_3 (km)')
set(gca,'FontSize',20,'LineWidth',1,'YDir','reverse')
%% use SVD to extract the important modes in postseismic timeseries
mpostvec = zeros(ss.N,length(lvec));
rpostvec = nan(length(lvec),1);

% extract time window
tind = t>0 & t<=700;% & t<=100;
ustamod = (usta(tind,:)-usta(find(tind,1),:) - ...
    repmat(Vpl/pi*atan(ox'./20e3),length(t(tind)),1).*repmat(t(tind),1,length(ox)));% remove deep creep signal
[U,S,V] = svd(ustamod,'econ');
Svals = diag(S);

% say we only want to use the sth modes
srank = 1;
% timeseries for the sth mode
% U(:,srank);
% spatial pattern
upost = V(:,srank);

% timeseries due to sth mode
uapprox = U(:,1:srank)*S(1:srank,1:srank)*V(:,1:srank)';

% invert for upost
for i = 1:length(lvec)
%     mpostvec(:,i) = inv(Gd'*Gd + lvec(i)^2*(Lsm'*Lsm))*Gd'*upost;
    mpostvec(:,i) = lsqnonneg([Gd;lvec(i).*Lsm],[upost;zeros(ss.N,1)]);
    rpostvec(i) = (upost-Gd*mpostvec(:,i))'*(upost-Gd*mpostvec(:,i));
end

figure(10)
subplot(313)
semilogx(lvec,rpostvec,'v-','LineWidth',2)
axis tight, grid on
ylabel('\epsilon^t \epsilon'), xlabel('\lambda')
set(gca,'FontSize',20,'LineWidth',1)

mst = 1; mend = 5;
% plot postseismic slip results
figure(21),clf
subplot(2,1,1)
plot(ox./1e3,upost,'ro'), hold on
plot(ox./1e3,Gd*mpostvec(:,mst:mend),'-','LineWidth',1)
axis tight, grid on
xlabel('x_2 (km)'), ylabel('Power of mode')
set(gca,'FontSize',20,'LineWidth',1)

subplot(2,1,2)
% plot(meqvec(:,mst:mend),ss.y3c./1e3,'.','LineWidth',1), hold on
plot(mpostvec(:,mst:mend),ss.y3c./1e3,'-','LineWidth',1)
axis tight, grid on
xlabel('Slip (mode)'), ylabel('x_3 (km)')
set(gca,'FontSize',20,'LineWidth',1,'YDir','reverse')
    
%% recombine to get slip distribution timeseries
mbest = 1;
spostapprox = Svals(srank).*mpostvec(:,mbest)*U(:,srank)';

figure(30),clf
subplot(2,1,1)
plot(ox./1e3,ueq,'ro','MarkerFaceColor',rgb('orange'),'MarkerEdgeColor',rgb('orange')), hold on
% plot(ox./1e3,ustamod(end,:),'bo')
plot(ox./1e3,uapprox(end,:),'d','MarkerFaceColor',rgb('forestgreen'),'MarkerEdgeColor',rgb('forestgreen'))

plot(ox./1e3,Gd*meqvec(:,mbest),'k-','LineWidth',1)
plot(ox./1e3,Gd*spostapprox(:,end),'k-','LineWidth',1)
axis tight, grid on
xlabel('x_2 (km)'), ylabel('u_1 (m)')
set(gca,'FontSize',20,'LineWidth',1)

subplot(2,1,2)
plot(meqvec(:,mbest),ss.y3c./1e3,'-','LineWidth',3,'Color',rgb('orange')), hold on
plot(spostapprox(:,end),ss.y3c./1e3,'-','LineWidth',2,'Color',rgb('forestgreen'))
plot(mintervec(:,mbest)./Vpl,ss.y3c./1e3,'-','LineWidth',2,'Color',rgb('steelblue'))
axis tight, grid on
% legend('coseismic','afterslip')
xlabel('Slip (m)'), ylabel('x_3 (km)')
set(gca,'FontSize',20,'LineWidth',1,'YDir','reverse')
    
%% derive stress-displacement relationship    
vinit = spostapprox(:,2) - spostapprox(:,1);
tvec = t(tind);

% fault stress
taupost = K12*meqvec(:,mbest) + K12*(spostapprox - Vpl*repmat(tvec',ss.N,1));

% compute v and deltau
vpost = abs(diff(spostapprox,1,2));
% dvpost = diff(spostapprox,2,2);
dvpost = diff(log(vpost),1,2);

dtaupost = diff((taupost),1,2);

figure(40),clf
plot(K12*meqvec(:,mbest),ss.y3c./1e3,'b-','LineWidth',2), hold on
plot(1000*vinit,ss.y3c./1e3,'r-','LineWidth',2)
axis tight, grid on
set(gca,'FontSize',20,'LineWidth',1,'YDir','reverse')

% estimate frictional properties
ind = [71:2:100];
figure(41),clf
for i = 1:length(ind)
%     subplot(211)
%     scatter(dvpost(ind(i),:),dtaupost(ind(i),2:end),100,tvec(3:end),'o','filled'), hold on
%     %xplt = linspace(min(dvpost(ind,:)),max(dvpost(ind,:)))';
%     % plot(xplt,xplt.*-.55/100,'r-','LineWidth',2)
%     % semilogx(dvpost(1,:),vpost(1,2:end)*5.4,'ro')
%     axis tight, grid on, box on
%     ylabel('\Delta\tau (MPa)'), xlabel('\Delta log(v)')
%     cb=colorbar;cb.Label.String = 'Time (days)';
%     set(gca,'FontSize',20,'LineWidth',1)
%     
%     subplot(212)
    toplot = (taupost(ind(i),2:end) - (min(taupost(ind(i),:))+.5*i)*1);
    scatter(vpost(ind(i),:),toplot,60,tvec(2:end),'o','filled','MarkerEdgeColor','none'), hold on
    axis tight, grid on, box on
    ylabel('\tau (MPa)'), xlabel('v (m/day)')
    cb=colorbar;cb.Label.String = 'Time (days)';
    set(gca,'FontSize',20,'LineWidth',1,'Xscale','log')%,'Yscale','lin')
end
    
    





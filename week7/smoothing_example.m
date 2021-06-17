% smoothing/regularization and hyperparameters in inverse problems
% Rishav Mallick

clear
addpath ~/Dropbox/scripts/utils

eq_dat = readtable('eq_displacement.dat');

% extract important variables
ox = eq_dat.ox;
ueq = eq_dat.u1_eq + 0.0005*randn(length(ox),1);% add noise for abic

 
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
%% run inversion with smoothing

lvec = 10.^[-6:.05:0]';
meqvec = zeros(ss.N,length(lvec));
resvec = zeros(length(lvec),1);
abicvec = zeros(length(lvec),1);

% define smoothing matrix
Lsm = compute_laplacian1d(ss.N);

for i = 1:length(lvec)
%     meqvec(:,i) = inv(Gd'*Gd + lvec(i)^2*(Lsm'*Lsm))*Gd'*ueq;
    meqvec(:,i) = lsqnonneg([Gd;lvec(i).*Lsm],[ueq;zeros(ss.N,1)]);
    resvec(i) = (ueq-Gd*meqvec(:,i))'*(ueq-Gd*meqvec(:,i));
    abicvec(i) = abic_alphabeta_sc(ueq,1,Gd,lvec(i),Lsm,0);
end

% choose models that are acceptable
loptind = find(abicvec==min(abicvec));
mst = loptind-2; mend = loptind+3;
cspec = jet(mend-mst+1);

% smoothing results/stats
figure(10),clf
subplot(211)
semilogx(lvec,resvec,'-','LineWidth',2), hold on
semilogx(lvec(mst:mend),resvec(mst:mend),'rp','MarkerSize',10)
axis tight, grid on
ylabel('\epsilon^t \epsilon'), xlabel('\lambda')
set(gca,'FontSize',20,'LineWidth',1)
  
subplot(212)
semilogx(lvec,abicvec-min(abicvec(:)),'-','LineWidth',2)
axis tight, grid on
ylabel('\Delta ABIC'), xlabel('\lambda')
set(gca,'FontSize',20,'LineWidth',1)

% plot slip results
figure(11),clf
subplot(2,1,1)
plot(ox./1e3,ueq,'ro'), hold on
for i = mst:mend
    plot(ox./1e3,Gd*meqvec(:,i),'-','LineWidth',1,'Color',cspec(i-mst+1,:))
end
axis tight, grid on
xlabel('x_2 (km)'), ylabel('u_1 (m)')
set(gca,'FontSize',20,'LineWidth',1)

subplot(2,1,2)
for i = mst:mend
    plot(meqvec(:,i),ss.y3c./1e3,'-','LineWidth',2,'Color',cspec(i-mst+1,:)), hold on
end
axis tight, grid on
xlabel('Slip (m)'), ylabel('x_3 (km)')
set(gca,'FontSize',20,'LineWidth',1,'YDir','reverse')

%% use Tikhonov 0th order

% define regularization matrix
Lsm = eye(ss.N);

for i = 1:length(lvec)
    %meqvec(:,i) = inv(Gd'*Gd + lvec(i)^2*(Lsm'*Lsm))*Gd'*ueq;
    meqvec(:,i) = lsqnonneg([Gd;lvec(i).*Lsm],[ueq;zeros(ss.N,1)]);
    resvec(i) = (ueq-Gd*meqvec(:,i))'*(ueq-Gd*meqvec(:,i));
    abicvec(i) = abic_alphabeta_sc(ueq,1,Gd,0,Lsm,lvec(i));
end

% choose models that are acceptable
loptind = find(abicvec==min(abicvec));
mst = loptind-2; mend = loptind+3;
cspec = jet(mend-mst+1);

figure(20),clf
subplot(211)
semilogx(lvec,resvec,'-','LineWidth',2), hold on
semilogx(lvec(mst:mend),resvec(mst:mend),'rp','MarkerSize',10)
axis tight, grid on
ylabel('\epsilon^t \epsilon'), xlabel('\beta')
set(gca,'FontSize',20,'LineWidth',1)
  
subplot(212)
semilogx(lvec,abicvec-min(abicvec(:)),'-','LineWidth',2), hold on
axis tight, grid on
ylabel('\Delta ABIC'), xlabel('\beta')
set(gca,'FontSize',20,'LineWidth',1)

% plot slip results
figure(21),clf
subplot(2,1,1)
plot(ox./1e3,ueq,'ro'), hold on
for i = mst:mend
    plot(ox./1e3,Gd*meqvec(:,i),'-','LineWidth',1,'Color',cspec(i-mst+1,:))
end
axis tight, grid on
xlabel('x_2 (km)'), ylabel('u_1 (m)')
set(gca,'FontSize',20,'LineWidth',1)

subplot(2,1,2)
for i = mst:mend
    plot(meqvec(:,i),ss.y3c./1e3,'-','LineWidth',2,'Color',cspec(i-mst+1,:)), hold on
end
axis tight, grid on
xlabel('Slip (m)'), ylabel('x_3 (km)')
set(gca,'FontSize',20,'LineWidth',1,'YDir','reverse')

%% use ABIC for both smoothing methods combined
Lsm = compute_laplacian1d(ss.N);

lvec = [logspace(-5,-3.5,5),logspace(-3.4,-.5,50)];
bvec = [logspace(-5,-4,5),logspace(-3.9,-2.5,40),logspace(-2.4,-2,5)];

[Lv,Bv] = meshgrid(lvec,bvec);
abicvec = zeros(size(Lv));
meqvec = zeros(ss.N,length(Lv(:)));

for i = 1:length(abicvec(:))
    %meqvec(:,i) = (Gd'*Gd + Lv(i)^2*(Lsm'*Lsm) + Bv(i)^2*eye(ss.N))\Gd'*ueq;
    meqvec(:,i) = lsqnonneg([Gd;Lv(i).*Lsm;Bv(i).*eye(ss.N)],[ueq;zeros(2*ss.N,1)]);
    abicvec(i) = abic_alphabeta(ueq,1,Gd,Lv(i),Lsm,Bv(i));
end

% choose models in a small range about the optimum model
moptind = find(abicvec==min(abicvec(:)));
modelrangeind = find(abicvec-min(abicvec(:))<=1);

figure(30),clf
pcolor(lvec,bvec,abicvec-min(abicvec(:))), hold on, %shading flat
contour(lvec,bvec,abicvec-min(abicvec(:)),[0:1:10],'w-','LineWidth',1)
axis tight, box on
cb=colorbar;cb.Label.String = '\Delta ABIC';cb.TickDirection='out';cb.LineWidth=2;
caxis([0 40])
colormap hot(32)
xlabel('\lambda'), ylabel('\beta')
set(gca,'YScale','log','XScale','log','FontSize',20,'LineWidth',2,'TickDir','both')

% plot results
cspec = jet(length(modelrangeind));
figure(31),clf
subplot(2,1,1)
plot(ox./1e3,ueq,'ro'), hold on
for i = 1:length(modelrangeind)
    plot(ox./1e3,Gd*meqvec(:,modelrangeind(i)),'-','LineWidth',1,'Color',cspec(i,:))
end
axis tight, grid on
xlabel('x_2 (km)'), ylabel('u_1 (m)')
set(gca,'FontSize',20,'LineWidth',1)

subplot(2,1,2)
for i = 1:length(modelrangeind)
    plot(meqvec(:,modelrangeind(i)),ss.y3c./1e3,'-','LineWidth',1,'Color',cspec(i,:)), hold on
end
axis tight, grid on
xlabel('Slip (m)'), ylabel('x_3 (km)')
set(gca,'FontSize',20,'LineWidth',1,'YDir','reverse')

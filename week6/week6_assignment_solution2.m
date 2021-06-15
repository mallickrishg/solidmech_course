% solutions to week 6 assignment part 2
% Rishav Mallick, EOS, 2021

clear
addpath ../unicycle/matlab/
import unicycle.*

%% load fault geometry and data

patchfname = 'flt_2D.seg';
earthModel = unicycle.greens.okada92(30e3,1/4);
rcv = unicycle.geometry.receiver(patchfname,earthModel);

% down-dip distance
dipdist = cumsum(rcv.W);

% load data
dat_eq = readtable('week6assignment/eq_planestrain.dat');
du = [dat_eq.ue;dat_eq.uz];

% compute displacement kernels
oxpred = linspace(-201e3,200e3,1000)';
[~,Gdip] = rcv.displacementKernels([dat_eq.ox,dat_eq.ox.*0,dat_eq.ox.*0],3);
[Gstrikepred,Gdippred] = rcv.displacementKernels([oxpred,oxpred.*0,oxpred.*0],3);
Ge = Gdip(1:3:end,:);
Gz = Gdip(3:3:end,:);

Gep = Gdippred(1:3:end,:);
Gzp = Gdippred(3:3:end,:);
Gnp = Gstrikepred(2:3:end,:);
Gmod = [Ge;Gz];

% linear inversion
% m1 = Gmod\du;
m1 = lsqnonneg(Gmod,du);

upred = Gmod*m1;
uepred = upred(1:end/2);
uzpred = upred(end/2+1:end);

figure(1),clf
plot(dat_eq.ox./1e3,dat_eq.ue,'s','MarkerFaceColor',rgb('forestgreen'),'MarkerEdgeColor','none'), hold on
plot(dat_eq.ox./1e3,dat_eq.uz,'s','MarkerFaceColor',rgb('royalblue'),'MarkerEdgeColor','none')
axis tight, grid on, box on
plot(oxpred./1e3,Gep*m1,'r-','LineWidth',2)
plot(oxpred./1e3,Gzp*m1,'r-','LineWidth',2)
xlabel('x_2 (km)'), ylabel('u (m)')
set(gca,'FontSize',20,'LineWidth',2)

figure(2),clf
plot(rcv.xc(:,1)./1e3,m1,'rx-','LineWidth',2), hold on
plot(oxpred,0.*oxpred,'k-','LineWidth',1)
axis tight, grid on
xlim([min(oxpred) max(oxpred)]./1e3)
xlabel('x_2 (km)'), ylabel('Slip (m)')
set(gca,'FontSize',20,'LineWidth',2)

% compute traction kernel
[~,~,~,Kdd,~,~] = rcv.tractionKernels(rcv);
% stress change on fault
tauchange = Kdd*m1;

figure(3),clf
subplot(211)
plot(dipdist./1e3,m1,'rx-','LineWidth',2)
axis tight, box on
xlabel('\zeta (km)'), ylabel('Slip (m)')
set(gca,'FontSize',20,'LineWidth',2)

subplot(212)
plot(dipdist./1e3,tauchange,'rx-','LineWidth',2), hold on
plot(dipdist./1e3,0.*dipdist,'k-','LineWidth',1)
axis tight, box on
xlabel('\zeta (km)'), ylabel('\Delta\tau (MPa)')
set(gca,'FontSize',20,'LineWidth',2)

figure(10),clf
mdummy = zeros(rcv.N,1);
mdummy(dipdist>10e3 & dipdist<50e3) = 4;
% mdummy=  m1;
plot(oxpred./1e3,Gnp*mdummy,'r-','LineWidth',2)
axis tight, grid on, box on
xlabel('x_2 (km)'), ylabel('u_n (m)')
%% interseismic velocity field
mparams = [1,rcv.dip(1),20e3,100e3];

[ueint,uzint] = espm2d(mparams,oxpred);

figure(4),clf
plot(oxpred./1e3,ueint,'r-','LineWidth',2), hold on
plot(oxpred./1e3,uzint,'-','LineWidth',2,'Color',rgb('steelblue'))
axis tight, box on, grid on
xlabel('x_2 (km)'), ylabel('u (m)')
set(gca,'FontSize',20,'LineWidth',2)



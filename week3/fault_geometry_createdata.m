% script to generate fault slip and traction data for class exercise
% Rishav Mallick, EOS, 2021

clear

addpath ~/Dropbox/scripts/unicycle/matlab/
import unicycle.*

%% read geometry file
patchfname = 'flt_3D.seg';
earthModel = unicycle.greens.okada92(30e3,1/4);
rcv = unicycle.geometry.receiver(patchfname,earthModel);

% compute traction kernels
[Kss,Kds,Ksd,Kdd,~,~]=rcv.tractionKernels(rcv);

%% compute slip 
rng(42)
tau_s = 0.5.*ones(rcv.N,1);slip_s = zeros(rcv.N,1);
tau_d = ones(rcv.N,1);slip_d = zeros(rcv.N,1);

pinfaultvec = randsample(rcv.N,40);

pinfault = true(rcv.N,1);
pinfault(pinfaultvec) = false;
pinfault(rcv.x(:,3)==0|rcv.xc(:,3) == min(rcv.xc(:,3))) = false;

tauvec = [tau_s(pinfault);tau_d(pinfault)];
Kmod = [Kss(pinfault,pinfault),Kds(pinfault,pinfault);Ksd(pinfault,pinfault),Kdd(pinfault,pinfault)];

slip = -Kmod\tauvec;
slip_s(pinfault) = slip(1:end/2);
slip_d(pinfault) = slip(end/2+1:end);

% calculate actual stress on fault
tau_fs = [Kss Kds]*[slip_s;slip_d];
tau_fd = [Ksd Kdd]*[slip_s;slip_d];

% plot
figure(1),clf
rcv.plotPatch(sqrt(slip_s.^2 + slip_d.^2)), hold on
rcv.plotSlipVectors(slip_s,slip_d,7e2,'b');
axis tight equal, box on, grid on
cb=colorbar;cb.Location = 'northoutside';cb.Label.String = 'Slip (m)';
view(12,28)
colormap(hot(100))
caxis([0 3])
xlabel('East (m)'), ylabel('North (m)'), zlabel('Depth (m)')
set(gca,'FontSize',20)

% calculate work done and average stress drop
dA = rcv.L.*rcv.W;

work = abs((tau_fs'*(slip_s.*dA) + tau_fd'*(slip_d.*dA))); % in MJ
tau_av = (tau_fs'*(slip_s.*dA) + tau_fd'*(slip_d.*dA))/sum(sqrt(slip_s.^2 + slip_d.^2).*dA);% in MPa

T = table(rcv.xc,rcv.L,rcv.W,rcv.strike,rcv.dip,slip_s,slip_d,tau_fs,tau_fd,'VariableNames',...
    {'xyz_c','L','W','strike','dip','slip_strike','slip_dip','traction_strike','traction_dip'});

writetable(T,'coseismicslipstress.dat','WriteVariableNames',1,'Delimiter','\t')
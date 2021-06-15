% script to read fault data and estimate fault geometry struct
% compute work done and average stress drop 
% Rishav Mallick, EOS, 2021

clear

datatab = readtable('coseismicslipstress.dat');

% store data in a structure
flt.xc = datatab{:,1:3};
flt.L = datatab.L;
flt.W = datatab.W;
flt.strike = datatab.strike;
flt.dip = datatab.dip;

slip_s = datatab.slip_strike;
slip_d = datatab.slip_dip;
tau_s = datatab.traction_strike;
tau_d = datatab.traction_dip;

slipmag = sqrt(datatab.slip_dip.^2 + datatab.slip_strike.^2);

flt.N = length(flt.L);

%% unit vectors
flt.nv = [cosd(flt.strike).*sind(flt.dip),-sind(flt.strike).*sind(flt.dip),cosd(flt.dip)];
flt.sv = [sind(flt.strike),cosd(flt.strike),zeros(flt.N,1)];
flt.dv = [cosd(flt.strike).*cosd(flt.dip),-sind(flt.strike).*cosd(flt.dip),-sind(flt.dip)];

%% plotting
flt.xTL = flt.xc - 0.5.*repmat(flt.L,1,3).*flt.sv - 0.5.*repmat(flt.W,1,3).*flt.dv;
flt.xTR = flt.xc + 0.5.*repmat(flt.L,1,3).*flt.sv - 0.5.*repmat(flt.W,1,3).*flt.dv;
flt.xBL = flt.xc - 0.5.*repmat(flt.L,1,3).*flt.sv + 0.5.*repmat(flt.W,1,3).*flt.dv;
flt.xBR = flt.xc + 0.5.*repmat(flt.L,1,3).*flt.sv + 0.5.*repmat(flt.W,1,3).*flt.dv;


figure(1),clf
for i = 1:flt.N
    plot3([flt.xTL(i,1),flt.xTR(i,1),flt.xBR(i,1),flt.xBL(i,1),flt.xTL(i,1)]./1e3,[flt.xTL(i,2),flt.xTR(i,2),flt.xBR(i,2),flt.xBL(i,2),flt.xTL(i,2)]./1e3,...
        [flt.xTL(i,3),flt.xTR(i,3),flt.xBR(i,3),flt.xBL(i,3),flt.xTL(i,3)]./1e3,'k-','LineWidth',2), hold on
end
axis tight equal, grid on, box on
view(15,30)
xlabel('East (km)'), ylabel('North (km)'), zlabel('Depth (km)')
set(gca,'FontSize',20)

figure(2),clf
for i = 1:flt.N
    patch([flt.xTL(i,1),flt.xTR(i,1),flt.xBR(i,1),flt.xBL(i,1),flt.xTL(i,1)]./1e3,[flt.xTL(i,2),flt.xTR(i,2),flt.xBR(i,2),flt.xBL(i,2),flt.xTL(i,2)]./1e3,...
        [flt.xTL(i,3),flt.xTR(i,3),flt.xBR(i,3),flt.xBL(i,3),flt.xTL(i,3)]./1e3,slipmag(i)), hold on
end

scf = 1e0;
slipvec = flt.sv.*repmat(datatab.slip_strike,1,3) - flt.dv.*repmat(datatab.slip_dip,1,3);
rake = atan2d(slip_d,slip_s);

% plot slipvectors
quiver3(flt.xc(:,1)./1e3,flt.xc(:,2)./1e3,flt.xc(:,3)./1e3,slipvec(:,1).*scf,slipvec(:,2).*scf,slipvec(:,3).*scf,'b-','LineWidth',2)

cb=colorbar;cb.Label.String = 'Slip (m)';
axis tight equal, grid on, box on
view(25,35)
caxis([0 3])
colormap(hot(100))
xlabel('East (km)'), ylabel('North (km)'), zlabel('Depth (km)')
set(gca,'FontSize',20)

%% plot tractions
% calculate work done and average stress drop
dA = flt.L.*flt.W;

work = abs((tau_s'*(slip_s.*dA) + tau_d'*(slip_d.*dA))); % in MJ
tau_av = (tau_s'*(slip_s.*dA) + tau_d'*(slip_d.*dA))/sum(sqrt(slip_s.^2 + slip_d.^2).*dA);% in MPa

% calculate stress change in the rake direction
tau_r = sum([tau_s tau_d].*[cosd(rake) sind(rake)],2);

figure(3),clf
plot_faultprop(flt,tau_r)
axis tight equal, grid on, box on
view(25,35)
cb=colorbar;cb.Label.String = '\tau_{rake} (MPa)';
caxis([-1 1].*abs(tau_av))
colormap bluewhitered
xlabel('East (km)'), ylabel('North (km)'), zlabel('Depth (km)')
set(gca,'FontSize',20)





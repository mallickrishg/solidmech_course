% script to compare displacements from fault slip
% comparing different tapering length-scales
% Rishav Mallick, EOS, 2021

clear

addpath ~/Dropbox/scripts/unicycle/matlab/
import unicycle.*

%% read geometry file
patchfname = 'flt_3D.seg';
earthModel = unicycle.greens.okada92(30e3,1/4);
rcv = unicycle.geometry.receiver(patchfname,earthModel);

% down-dip distance
dipdist = cumsum(rcv.W);

nobs = 1000;
ox = linspace(-20e3,100e3,nobs)';
[~,Gd] = rcv.displacementKernels([ox ox.*0 ox.*0],3);

Gde = Gd(1:3:end,:);
Gdz = Gd(3:3:end,:);

rcv.slip = ones(rcv.N,1);

% pin some regions
xpin = [10e3, 80e3];
% linear taper width
tw = linspace(0e3,10e3,20);
cspec = jet(length(tw));
% ind = rcv.xc(:,1)<=xpin(1) | rcv.xc(:,1)>=xpin(2);

figure(1),clf
figure(2),clf
for i = 1:length(tw)
    ind = dipdist<=xpin(1) | dipdist>=xpin(2);
    rcv.slip(ind) = 0;
    % add shallow linear taper
    ind = dipdist<=(xpin(1)+tw(i)) & dipdist>=xpin(1);
    ntaper = length(find(ind));
    rcv.slip(ind) = linspace(0,1,ntaper);
    
    % calculated displacements
    uz = Gdz*rcv.slip;

    figure(1)
    subplot(211)
    plot(ox./1e3,uz,'-','LineWidth',2,'Color',cspec(i,:)), hold on    
    axis tight, grid on
    set(gca,'Fontsize',15)
    ylim([-0.1 0.6])
    xlabel('x_2 (km)'), ylabel('u_z (m)')
    
    subplot(212)
    plot(rcv.xc(:,1)./1e3,rcv.slip,'-','LineWidth',2,'Color',cspec(i,:)), hold on
    axis tight, grid on
    set(gca,'Fontsize',15)
    ylim([-0.1 1.1])
    xlim([min(ox./1e3),max(ox./1e3)])
    xlabel('x_2 (km)'), ylabel('Slip (m)')
    
    figure(2)
    plot(rcv.slip,rcv.xc(:,3)./1e3,'-','LineWidth',2,'Color',cspec(i,:)), hold on
    axis tight, grid on
    set(gca,'Fontsize',15)
    ylabel('Depth (km)'), xlabel('Slip (m)')
%     yyaxis right
%     plot(rcv.slip,-rcv.xc(:,1)./1e3,'--','LineWidth',1,'Color',cspec(i,:)), hold on
%     axis tight, grid on
%     set(gca,'Fontsize',15)
%     ylabel('Depth (km)'), xlabel('Slip (m)')
end
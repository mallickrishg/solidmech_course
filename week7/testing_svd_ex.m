% playing around with svd
% Rishav Mallick, EOS, 2021

clear
addpath ~/Dropbox/scripts/utils

eq_dat = readtable('assignment_week7/eq_displacement.dat');
gps = readtable('assignment_week7/gps_timeseries.dat');

% extract important variables
ox = eq_dat.ox;
t = gps.tvals;
usta = gps{:,2:end};

%% use svd and extract modes
tind = t>0 & t<20;
ustamod = usta(tind,:);

[U,S,V] = svd(ustamod,'econ');
Svals = diag(S);
srank = [1,2,10];

figure(1),clf
subplot(211)
semilogy(Svals,'k-','LineWidth',2)
axis tight, grid on

subplot(212)
plot(cumsum(Svals)./sum(Svals),'k-','LineWidth',2)
axis tight

figure(2),clf
subplot(2,2,1)
imagesc(ox./1e3,t(tind),ustamod), hold on
contour(ox./1e3,t(tind),ustamod,[-1:0.1:1].*0.1,'k-')
xlabel('x_2 (km)'), ylabel('t (days)')
cb=colorbar;
caxis([-1 1].*0.15)
colormap bluewhitered(100)
set(gca,'FontSize',15,'YDir','normal')

for i=1:length(srank)
    subplot(2,2,i+1)
    uapprox = U(:,1:srank(i))*S(1:srank(i),1:srank(i))*V(:,1:srank(i))';
    imagesc(ox./1e3,t(tind),uapprox), hold on
    contour(ox./1e3,t(tind),uapprox,[-1:0.1:1].*0.1,'k-')
    xlabel('x_2 (km)'), ylabel('t (days)')
    cb=colorbar;
    caxis([-1 1].*0.15)
    set(gca,'FontSize',15,'YDir','normal')
end
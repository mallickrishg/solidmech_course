% playing around with svd
% Rishav Mallick, EOS, 2021

clear
addpath ~/Dropbox/scripts/utils

eq_dat = readtable('eq_displacement.dat');
gps = readtable('gps_timeseries.dat');

% extract important variables
ox = eq_dat.ox;
t = gps.tvals;
usta = gps{:,2:end};

%% use svd and extract modes
% tind = t<0 & t>=-20;
tind = t>0 & t<=100;
ustamod = (usta(tind,:)-usta(find(tind,1),:))*1000;
cval = 4;

[U,S,V] = svd(ustamod,'econ');
Svals = diag(S);
srank = [1,2,5];

figure(10),clf
xind = ox>0;
xindvec = find(xind);
cspec = jet(length(xindvec));
subplot(2,2,1)
for j = 1:length(xindvec)
    plot(t(tind),ustamod(:,xindvec(j)),'-o','LineWidth',1,'Color',cspec(j,:)), hold on
end
axis tight, grid on
ylabel('u_1 (mm)'), xlabel('Time (days)')
set(gca,'FontSize',15)


figure(1),clf
subplot(211)
semilogy(Svals,'k-','LineWidth',2)
axis tight, grid on

subplot(212)
plot(cumsum(Svals)./sum(Svals),'k-','LineWidth',2)
axis tight

figure(2),clf
subplot(2,2,1)
imagesc(ox(xind)./1e3,t(tind),ustamod(:,xind)), hold on
contour(ox(xind)./1e3,t(tind),ustamod(:,xind),[-1:0.1:1].*cval,'k-')
xlabel('x_2 (km)'), ylabel('t (days)')
cb=colorbar;
caxis([0 1].*cval)
colormap parula(100)
set(gca,'FontSize',15,'YDir','normal')

for i=1:length(srank)
    figure(2)
    subplot(2,2,i+1)
    uapprox = U(:,1:srank(i))*S(1:srank(i),1:srank(i))*V(:,1:srank(i))';
    imagesc(ox(xind)./1e3,t(tind),uapprox(:,xind)), hold on
    contour(ox(xind)./1e3,t(tind),uapprox(:,xind),[-1:0.1:1].*cval,'k-')
    xlabel('x_2 (km)'), ylabel('t (days)')
    cb=colorbar;
    caxis([0 1].*cval)
    title(['Number of modes = ' num2str(srank(i))])
    set(gca,'FontSize',15,'YDir','normal')
    
    figure(10)
    subplot(2,2,i+1)
    for j = 1:length(xindvec)
        plot(t(tind),ustamod(:,xindvec(j)),'-o','LineWidth',1,'Color',cspec(j,:)), hold on
    end
    axis tight, grid on
    ylabel('u_1 (mm)'), xlabel('Time (days)')
    title(['Number of modes = ' num2str(srank(i))])
    set(gca,'FontSize',15)
end

%% plot singular values (spatial modes and temporal modes)
cspec = [rgb('steelblue');rgb('orange');rgb('peru');rgb('royalblue');rgb('black')];%parula(srank(end));
figure(3),clf

subplot(211)
for i = 1:srank(end)
    plot(ox./1e3,Svals(i).*V(:,i),'-o','LineWidth',2,'Color',cspec(i,:)), hold on
end
ylabel('Power'), xlabel('x_2 (km)')
axis tight, grid on
set(gca,'FontSize',15,'LineWidth',2)

subplot(212)
for i = 1:srank(end)
    plot(t(tind),Svals(i).*U(:,i),'-o','LineWidth',2,'Color',cspec(i,:)), hold on
end
ylabel('Power'), xlabel('t (days)')
axis tight, grid on
set(gca,'FontSize',15,'LineWidth',2)
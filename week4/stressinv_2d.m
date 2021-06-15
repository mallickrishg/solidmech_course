% simple script to compute 2-d deviatoric stress state from a number of
% slip vector observations using least squares
% Rishav Mallick, 2021, EOS

clear
rng(42)

nobs = 100;
dip = 35 + 30.*rand(nobs,1);

% data = readtable('slipvec2d.dat');
% dip = data{:,1};
% nobs = length(data{:,1});


% define nhat and dhat functions
n = @(dip) [sind(dip);cosd(dip)];

d = @(dip) [cosd(dip);-sind(dip)];
r = @(dip) [cosd(dip);-sind(dip)];% thrust sense


% design matrix related stress tensor to fault slip vector
Gfun = @(n) [n(1)-n(1)^3+n(1)*n(2)^2, n(2)-2*n(1)^2*n(2);...
    n(2)^3-n(2)-n(1)^2*n(2), n(1)-2*n(1)*n(2)^2];

G = zeros(nobs*2,2);
dhatobs = zeros(nobs*2,1);
rhatobs = zeros(nobs*2,1);
nhatobs = zeros(nobs*2,1);
for i = 1:nobs
    G(2*i-1:2*i,:) = Gfun(n(dip(i)));
    dhatobs(2*i-1:2*i,1) = d(dip(i));
    rhatobs(2*i-1:2*i,1) = r(dip(i));
    nhatobs(2*i-1:2*i,1) = n(dip(i));
end

% given fault slip observations, infer the stress direction
spred = G\rhatobs;
Cs = inv(G'*G);


% sample and create 2-d distribution
nsample = 1e4;
spsample = mvnrnd(spred,Cs,nsample);
thetasample = .5*atan2(spsample(:,2),spsample(:,1));

figure(10),clf
subplot(2,2,1)
histogram(spsample(:,2),'Normalization','pdf','EdgeColor','none');
axis tight, box on
camroll(90)
ylabel('p'),xlabel('\sigma_{23}')
set(gca,'FontSize',20)

subplot(2,2,4)
histogram(spsample(:,1),'Normalization','pdf','EdgeColor','none')
axis tight, box on
ylabel('p'),xlabel('\sigma_{22}')
set(gca,'FontSize',20)

subplot(2,2,2)
histogram2(spsample(:,1),spsample(:,2),'Normalization','pdf','DisplayStyle','tile','EdgeColor','none','ShowEmptyBins','on'), hold on
[hvals,XBinEdges,YBinEdges] = histcounts2(spsample(:,1),spsample(:,2),'Normalization','probability');
contour(XBinEdges(1:end-1),YBinEdges(1:end-1),hvals',10,'Color',rgb('white'),'Linewidth',1)
colormap hot
plot(spred(1),spred(2),'rp','MarkerFaceColor','r','MarkerSize',20)
axis tight, grid on, box on
xlabel('\sigma_{22}'),ylabel('\sigma_{23}')
set(gca,'FontSize',20)

subplot(2,2,3)
polarhistogram(thetasample,'Normalization','pdf','EdgeColor','none')
axis tight, box on
title('\theta (\sigma_{max})')
set(gca,'FontSize',20)


Evec1 = zeros(2,nsample);Evec2 = zeros(2,nsample);
Eval = zeros(2,nsample);
for i = 1:nsample
    s22p = spsample(i,1);
    s23p = spsample(i,2);

    % compute eigenvectors of the stress tensor
    [evec,eval] = eig([s22p,s23p;s23p,-s22p]);
    Evec1(:,i) = evec(:,1);
    Evec2(:,i) = evec(:,2);
    Eval(:,i) = diag(eval);
end

% % if you want to predict the slip vector from the stress state
% dpred = G*spred;
% dpred = dpred./sqrt(dpred'*dpred);

%% plot directions
figure(1),clf
% plot fault orientation (d and n)
for i = 1:nobs
    %plot([-1 1]*dhatobs(2*i-1),[-1 1]*dhatobs(2*i),'k-p','LineWidth',.1), hold on
    plot([0 1]*dhatobs(2*i-1),[0 1]*dhatobs(2*i),'k-p','LineWidth',.1), hold on
    %plot([-1 1]*nhatobs(2*i-1),[-1 1]*nhatobs(2*i),'k-x','MarkerFaceColor','k','LineWidth',1)
end

plot(Evec1(1,:),Evec1(2,:),'r.')
% plot(-Evec1(1,:),-Evec1(2,:),'ro')
plot(Evec2(1,:),Evec2(2,:),'b.')
plot(-Evec2(1,:),-Evec2(2,:),'b.')
plot([-1 1],[0 0],'k--','LineWidth',2)
axis tight equal, grid on, box on
% colormap(flipud(hot))
xlabel('x_2'),ylabel('x_3')
set(gca,'FontSize',20)

% calculation angles of s1,s2
theta1 = atan2(Evec1(2,:),Evec1(1,:))';
theta2 = atan2(Evec2(2,:),Evec2(1,:))';
%
figure(2),clf
polarhistogram(theta2,20,'Normalization','pdf','Edgecolor','none','Facealpha',.6), hold on
polarhistogram(theta1,20,'Normalization','pdf','Edgecolor','none','Facealpha',.6)
polarhistogram(deg2rad(-dip),'Normalization','pdf','DisplayStyle','stairs','Linewidth',3,'EdgeColor','k')
axis tight
legend('\sigma_2','\sigma_1','\delta')
set(legend,'box','off')
set(gca,'FontSize',20,'GridAlpha',.3,'LineWidth',1,'Box','on')





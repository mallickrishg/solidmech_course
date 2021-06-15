% script to compute 3-d deviatoric stress state from a number of
% slip vector observations using least squares
% Rishav Mallick, 2021, EOS

clear

n = @(strike,dip) [cosd(strike)*sind(dip);-sind(strike)*sind(dip);cosd(dip)];
% d = @(strike,dip) [cosd(strike)*cosd(dip);-sind(strike)*sind(dip);-sind(dip)];
r = @(strike,dip,rake) [cosd(rake)*sind(strike)-sind(rake)*cosd(strike)*cosd(dip);...
    cosd(rake)*cosd(strike)+sind(rake)*sind(strike)*cosd(dip);...
    sind(rake)*sind(dip)]; % original equation

% design matrix

Gfun = @(n) [n(1)-n(1)^3+n(1)*n(3)^2, n(2)-2*n(2)*n(1)^2, n(3)-2*n(3)*n(1)^2, -n(1)*n(2)^2+n(1)*n(3)^2, -2*n(1)*n(2)*n(3);...
            -n(2)*n(1)^2+n(2)*n(3)^2, n(1)-2*n(1)*n(2)^2, -2*n(1)*n(2)*n(3),   n(2)-n(2)^3+n(2)*n(3)^2,  n(3)-2*n(3)*n(2)^2;...
            -n(3)*n(1)^2-n(3)+n(3)^3, -2*n(1)*n(2)*n(3),  n(1)-2*n(1)*n(3)^2, -n(2)^2*n(3)-n(3)+n(3)^3,  n(2)-2*n(2)*n(3)^2];

%% input data
% strike = 45;
% dip = 60;
% rake = 0;
% % rinit = r(strike,dip,rake);
% nobs = 500;
% 
% strikeobs = strike + 10.*randn(nobs,1);
% dipobs = dip - 5 + 10.*rand(nobs,1);
% rakeobs = rake + 1.*randn(nobs,1);

% load('INPUT_example0D.mat');
% datamat = example0D;
% spredsatsi = [1.01345600000000;-0.506121000000000;-0.153006000000000;-0.0888350000000000;-0.279231000000000];
% spredsatsimat = [spredsatsi(1:3)';spredsatsi([2,4,5])';spredsatsi([3,5])',-(spredsatsi(1)+spredsatsi(4))];
% [evecsatsi,evalsatsi] = eig(spredsatsimat);

% load('INPUT_example1D.mat');
% ind = example1D(:,2)<=10;
% datamat = example1D(ind,:);

load('INPUTdata.mat');
datamat = T;
nobs = length(datamat(:,1));
strikeobs = datamat(:,3);
dipobs = datamat(:,4);
rakeobs = datamat(:,5);

%% create design matrix
G = zeros(nobs*3,5);
rhatobs = zeros(nobs*3,1);
nhatobs = zeros(nobs*3,1);

for i = 1:nobs
    G(3*i-2:3*i,:) = Gfun(n(strikeobs(i),dipobs(i)));
    rhatobs(3*i-2:3*i,1) = r(strikeobs(i),dipobs(i),rakeobs(i));
    nhatobs(3*i-2:3*i,1) = n(strikeobs(i),dipobs(i));
end


% order followed in stress tensor - [s11;s12;s13;s22;s23]
spred = G\rhatobs;
Cs = inv(G'*G);

spredmat = [spred(1:3)';spred([2,4,5])';spred([3,5])',-(spred(1)+spred(4))];
[evecopt,evalopt] = eig(spredmat);
evalopt = diag(evalopt);

% sample and create 2-d distribution
nsample = 1e4;
spsample = mvnrnd(spred,Cs,nsample);

Evec1 = zeros(3,nsample);
Evec2 = zeros(3,nsample);
Evec3 = zeros(3,nsample);
Eval = zeros(3,nsample);

for i = 1:nsample
    s11p = spsample(i,1);
    s12p = spsample(i,2);
    s13p = spsample(i,3);
    s22p = spsample(i,4);
    s23p = spsample(i,5);

    % compute eigenvectors of the stress tensor
    [evec,eval] = eig([s11p,s12p,s13p;s12p,s22p,s23p;s13p,s23p,-(s11p+s22p)]);
    Evec1(:,i) = evec(:,1);
    Evec2(:,i) = evec(:,2);
    Evec3(:,i) = evec(:,3);
    Eval(:,i) = diag(eval);
end

% calculate stress ratio
phi = ((Eval(2,:))-(Eval(3,:)))./((Eval(1,:))-(Eval(3,:)));
R = ((Eval(1,:))-(Eval(2,:)))./((Eval(1,:))-(Eval(3,:)));

% plot stress orientations as colored dots
figure(1),clf
for i = 1:nobs
    plot3(rhatobs(3*i-2),rhatobs(3*i-1),rhatobs(3*i),'kp'), hold on
    %plot3(nhatobs(3*i-2),nhatobs(3*i-1),nhatobs(3*i),'kx')
end
plot3(-sign(sign(Eval(1,:)).*Evec1(3,:)).*sign(Eval(1,:)).*Evec1(1,:),...
    -sign(sign(Eval(1,:)).*Evec1(3,:)).*sign(Eval(1,:)).*Evec1(2,:),...
    -sign(sign(Eval(1,:)).*Evec1(3,:)).*sign(Eval(1,:)).*Evec1(3,:),'.','Color',rgb('orange'))

plot3(-sign(sign(Eval(2,:)).*Evec2(3,:)).*sign(Eval(2,:)).*Evec2(1,:),...
    -sign(sign(Eval(2,:)).*Evec2(3,:)).*sign(Eval(2,:)).*Evec2(2,:),...
    -sign(sign(Eval(2,:)).*Evec2(3,:)).*sign(Eval(2,:)).*Evec2(3,:),'.','Color',rgb('forestgreen'))

plot3(-sign(sign(Eval(3,:)).*Evec3(3,:)).*sign(Eval(3,:)).*Evec3(1,:),...
    -sign(sign(Eval(3,:)).*Evec3(3,:)).*sign(Eval(3,:)).*Evec3(2,:),...
    -sign(sign(Eval(3,:)).*Evec3(3,:)).*sign(Eval(3,:)).*Evec3(3,:),'.','Color',rgb('skyblue'))

axis tight equal, box on, grid on
plot3(cosd(linspace(0,180)),sind(linspace(0,180)),zeros(100,1),'k-','LineWidth',2)
plot3(cosd(linspace(0,180)),-sind(linspace(0,180)),zeros(100,1),'k-','LineWidth',2)
plot3(zeros(100,1),sind(linspace(0,180)),cosd(linspace(0,180)),'k-','LineWidth',2)
plot3(zeros(100,1),-sind(linspace(0,180)),cosd(linspace(0,180)),'k-','LineWidth',2)
plot3(sind(linspace(0,180)),zeros(100,1),cosd(linspace(0,180)),'k-','LineWidth',2)
plot3(-sind(linspace(0,180)),zeros(100,1),cosd(linspace(0,180)),'k-','LineWidth',2)

% plot3(sign(evalopt').*evecopt(1,:),sign(evalopt').*evecopt(2,:),sign(evalopt').*evecopt(3,:),'k+','LineWidth',3,'MarkerSize',10)
% plot3(evecsatsi(1,:),evecsatsi(2,:),evecsatsi(3,:),'rx','LineWidth',3,'MarkerSize',10)

view([0 90])
xlabel('East'),ylabel('North')
set(gca,'FontSize',15)

%%%%%%%%% contours of posterior PDF
figure(2),clf
% plot rake vectors

[hvals,XBinEdges,YBinEdges] = histcounts2(-sign(sign(Eval(1,:)).*Evec1(3,:)).*sign(Eval(1,:)).*Evec1(1,:),...
    -sign(sign(Eval(1,:)).*Evec1(3,:)).*sign(Eval(1,:)).*Evec1(2,:),...
    [-1:0.05:1.05]',[-1:0.05:1.05]','Normalization','probability');
contour(XBinEdges(1:end-1),YBinEdges(1:end-1),hvals',10,'Color',rgb('orange'),'Linewidth',2), hold on            

[hvals,XBinEdges,YBinEdges] = histcounts2(-sign(sign(Eval(2,:)).*Evec2(3,:)).*sign(Eval(2,:)).*Evec2(1,:),...
    -sign(sign(Eval(2,:)).*Evec2(3,:)).*sign(Eval(2,:)).*Evec2(2,:),...
    [-1:0.05:1.05]',[-1:0.05:1.05]','Normalization','probability');
contour(XBinEdges(1:end-1),YBinEdges(1:end-1),hvals',10,'Color',rgb('forestgreen'),'Linewidth',1)            

[hvals,XBinEdges,YBinEdges] = histcounts2(-sign(sign(Eval(3,:)).*Evec3(3,:)).*sign(Eval(3,:)).*Evec3(1,:),...
    -sign(sign(Eval(3,:)).*Evec3(3,:)).*sign(Eval(3,:)).*Evec3(2,:),...
    [-1:0.05:1.05]',[-1:0.05:1.05]','Normalization','probability');
contour(XBinEdges(1:end-1),YBinEdges(1:end-1),hvals',10,'Color',rgb('skyblue'),'Linewidth',2)            

for i = 1:nobs
    plot(rhatobs(3*i-2),rhatobs(3*i-1),'kp','MarkerSize',3), hold on
end
plot(cosd(linspace(0,180)),sind(linspace(0,180)),'k-','LineWidth',2)
plot(cosd(linspace(0,180))./2,sind(linspace(0,180))./2,'k-','LineWidth',1)
plot(cosd(linspace(0,180)),-sind(linspace(0,180)),'k-','LineWidth',2)
plot(cosd(linspace(0,180))./2,-sind(linspace(0,180))./2,'k-','LineWidth',1)

thetavec = [0:30:360];
for i = 1:length(thetavec)
    plot([0 cosd(thetavec(i))],[0 sind(thetavec(i))],'k-','Linewidth',.1)
end

% plot(sign(evalopt').*evecopt(1,:),sign(evalopt').*evecopt(2,:),'k+','LineWidth',3,'MarkerSize',10)
% plot(evecsatsi(1,:),evecsatsi(2,:),'r+','LineWidth',3,'MarkerSize',10)

legend('\sigma_1','\sigma_2','\sigma_3','\lambda')
set(legend,'box','off')%,'location','northeast')
axis tight equal, box on, grid off
xlabel('East'),ylabel('North')
set(gca,'FontSize',20,'LineWidth',1,'YTick',[-1:0.5:1])

%% histogram of stress ratio
figure(3),clf
%plot phi
h1=histogram(phi,'Normalization','probability','EdgeColor','none');
hold on
hmax=h1.BinEdges(h1.Values==max(h1.Values))+h1.BinWidth/2;
plot([1 1]*hmax(1),get(gca,'ylim'),'b','lineWidth',2);

% plot R
h1=histogram(R,'Normalization','probability','EdgeColor','none');
hmax=h1.BinEdges(h1.Values==max(h1.Values))+h1.BinWidth/2;
plot([1 1]*hmax(1),get(gca,'ylim'),'r','lineWidth',2);

legend('\phi','\phi_{max}','R','R_{max}')
set(legend,'box','off','location','best')
xlabel('\phi or R'), ylabel('PDF')
axis tight, box on, grid on
xlim([0 1])
set(gca,'FontSize',15,'LineWidth',1,'YGrid','off')





















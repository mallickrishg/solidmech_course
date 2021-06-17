% AUTHOR:
% Rishav Mallick, 2019, Earth Observatory of Singapore

% Simulates earthquake cycles on 
% strike-slip faults in 2D antiplane.

clear
close all
addpath ~/Dropbox/scripts/utils/

% Rigidity (MPa)
G = 30e3; %MPa

% fault patches
M = 300; % patch size = 20e3/M

y2i = 0;
y3i = zeros(length(y2i),1);

% plate velocity
Vpl = 1e-9; %(m/s)
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                        M E S H                       %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % %  %
ss = [];
ss.G = G;
% Brittle-Ductile Tranisition Depth (m)
Transition = 20e3;
ss.y3f = zeros(length(y2i)*M,1);
ss.y2f = ss.y3f;
ss.Wf = ss.y3f;
ss.M = length(y2i)*M;

for i = 1:length(y2i)
    % Fault Meshes
    
    % top of fault patches
    y3f = linspace(y3i(i),Transition,M+1)';
    ss.y3f(1+(i-1)*M:i*M) = y3f(1:end-1);
    
    ss.y2f(1+(i-1)*M:i*M) = y2i(i).*ones(M,1);
    
    % width of fault patches
    ss.Wf(1+(i-1)*M:i*M)     = diff(y3f);
    
end

% mid of fault patches
ss.y3c = ss.y3f + ss.Wf./2;
ss.y2c = ss.y2f;

% compute traction kernels
evl = compute_stresskernels(ss);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         F R I C T I O N   P A R A M E T E R S        %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% effective confining pressure on fault (MPa)
ss.sigma = 50.*ones(ss.M,1);
% frictional parameters 
ss.a = 1e-2*ones(ss.M,1);
ss.b = ss.a - 5e-3*ones(ss.M,1);
% region of velocity-weakening
vw = ss.y3c>=2e3 & ss.y3c<=15e3;
ss.b(vw) = ss.a(vw) + 5e-3;
% static friction coefficient
ss.mu0 = 0.6*ones(ss.M,1);
% characteristic weakening distance (m)
ss.l = 0.004*ones(ss.M,1);
% plate velocity (m/s)
ss.Vpl = Vpl*ones(ss.M,1);
% reference slip rate (m/s)
ss.Vo = 1e-6*ones(ss.M,1);
% shear wave speed (m/s)
ss.Vs = 3e3*ones(ss.M,1);
% Degrees of Freedom [slip, log(v/vo), log(vo theta/L)]
ss.dgf = 3;
disp('Assigned Frictional Properties')
% Loading Stresses
evl.tau0 = -evl.K*ss.Vpl;
% discretization criteria
hstar_rr = pi/4.*G*ss.l./(ss.b-ss.a)./ss.sigma;
hstar_ra = pi/2.*G*ss.b.*ss.l./(ss.b-ss.a).^2/ss.sigma;
hstar_cz = hstar_rr.*(ss.b-ss.a)./ss.b;
hstar_rr = min(hstar_rr(vw));
hstar_ra = min(hstar_ra(vw));
hstar_cz = min(hstar_cz(vw));
dx_fault = (max(ss.Wf));
disp(['h*rr = ' num2str(hstar_rr) 'm'])
disp(['h*ra = ' num2str(hstar_ra) 'm'])
disp(['h*cz = ' num2str(hstar_cz) 'm'])
disp(['dx(fault) = ' num2str(dx_fault) 'm'])
disp(['Instability ratio L_vw/h*rr = ' num2str(sum(vw.*ss.Wf)/hstar_rr)])

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
%         N U M E R I C A L   S O L U T I O N           %
%                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
disp('Solving ODE')

% Initialize State Vector
Y0=zeros(ss.M*ss.dgf,1);
% Fault patches
Y0(1:ss.dgf:ss.M*ss.dgf) = zeros(ss.M,1);% slip
Y0(2:ss.dgf:ss.M*ss.dgf) = log(ss.Vo./ss.Vpl); % th
vinit = log(ss.Vpl*0.99./ss.Vo);
Y0(3:ss.dgf:ss.M*ss.dgf) = vinit; % v

% Simulation Duration
Tend = 20/Vpl/3.15e7;%years

tic
% initialize the function handle with
% set constitutive parameters
yp=@(t,y) odeElastic_a(t,y,ss,evl);
tic
% Solve the system
options=odeset('Refine',1,'RelTol',1e-6,'InitialStep',1e-6,'MaxStep',1e10); 
[t,Y]=ode45(yp,[0 Tend*3.15e7],Y0,options);
toc

f=repmat(ss.mu0',size(Y,1),1) + repmat(ss.a',size(Y,1),1).*Y(:,3:ss.dgf:ss.M*ss.dgf) + repmat(ss.b',size(Y,1),1).*Y(:,2:ss.dgf:ss.M*ss.dgf);
V=repmat(ss.Vo',size(Y,1),1).*exp(Y(1:end,3:ss.dgf:ss.M*ss.dgf));
slip = Y(:,1:ss.dgf:ss.M*ss.dgf);

%% extract earthquake events
[tvec_eq,~,slip_eq,L] = compute_earthquakestatistics(ss,logical(ones(ss.M,1)),t,V,slip);
Lsum = sum(L,2);

%% plot velocity timeseries
cspec=[cmap('steelblue',100,10,10);flipud(cmap('orange',50,57,10));flipud(cmap('orangered',50,20,25))];

xl = [1 Tend];
xlt = [find(abs(t./3.15e7-xl(1))==min(abs(t./3.15e7-xl(1)))) find(abs(t./3.15e7-xl(2))==min(abs(t./3.15e7-xl(2))))];
tdown = 5;

figure(1),clf
subplot(4,1,[1,2,3])
pcolor((t(xlt(1):tdown:xlt(2)) - t(xlt(1)))./3.15e7,ss.y3c(:)./1e3,(V(xlt(1):tdown:xlt(2),:))'./Vpl), shading interp, box on
hold on
for ii = 1:length(tvec_eq)
    pos1 = find(L(ii,:),1,'first');
    pos2 = find(L(ii,:),1,'last');
    plot(tvec_eq(ii)./3.15e7.*[1 1],ss.y3c([pos1 pos2])./1e3,'r-','LineWidth',1)
end
set(gca,'YDir','reverse','FontSize',15,'ColorScale','log')
xlabel('Time (yrs)')
colormap(cspec)
caxis(10.^[-2 2])
ylabel('\zeta_d (km)')

cb=colorbar;
cb.Location = 'northoutside';cb.Label.String = 'normalized V';

subplot(4,1,4)
index = ss.a<ss.b;
plot(t(xlt(1):xlt(2))./3.15e7,max(V(xlt(1):xlt(2),index),[],2),'-','LineWidth',2)
grid on, axis tight
ylim(10.^[-10 1])
xlabel('Time (yrs)'), ylabel('V (m/s)')
set(gca,'YScale','log','FontSize',15)

%% surface observations
% choose earthquake
neq = 12;

% resample data to mimic GPS
tresamp = (tvec_eq(neq)-3.15e7*2:86400:tvec_eq(neq)+3.15e7*2)';
Vresamp = interp1(t,V,tresamp);
slipresamp = interp1(t,slip,tresamp);

nobs = 60;
obs = [linspace(-100e3,100e3,nobs)' zeros(nobs,1)];

% displacement kernels
Gd = compute_displacementkernels(obs,ss);

% compute displacements and disp-rates
veq = (Gd.kd*Vresamp')';
ueq = (Gd.kd*slipresamp')';

v1 = veq + repmat(Vpl/pi*atan2(obs(:,1),Transition)',length(tresamp),1);
u1 = ueq + repmat(Vpl/pi*atan2(obs(:,1),Transition)',length(tresamp),1).*repmat(tresamp,1,nobs);

u1 = u1 - repmat(u1(1,:),length(tresamp),1);

% export data

u1_eq = Gd.kd*slip_eq(neq,:)';

figure(3),clf
subplot(131)
plot(slip_eq(neq,:),ss.y3c./1e3,'r-','LineWidth',3)
axis tight
xlabel('Slip (m)'), ylabel('Depth (km)')
set(gca,'Fontsize',15,'YDir','reverse','LineWidth',1)

subplot(1,3,[2,3])
plot(obs(:,1)./1e3,u1_eq,'rd','LineWidth',1,'MarkerFaceColor','r','MarkerSize',10)
axis tight, grid on
% ylim([0 max(u1_eq)]), xlim([0 max(obs(:,1)./1e3)])
ylabel('u_1 (m)'), xlabel('x_2 (km)')
set(gca,'Fontsize',15,'LineWidth',1)
%%
ox = obs(:,1);

T = table(ox,u1_eq);
% writetable(T,'assignment_week7/eq_displacement.dat');
writetable(T,'eq_displacement.dat');
%
figure(2),clf
subplot(211)
plot(tresamp./3.15e7, u1(:,1:4:end),'.','LineWidth',1)
axis tight, grid on
xlabel('Time (yr)'), ylabel('u_1 (m)')
set(gca,'Fontsize',15)

subplot(212)
plot(tresamp./3.15e7, v1(:,1:4:end)./Vpl,'.','LineWidth',1)
axis tight, grid on
ylim([0 1]*2)
xlabel('Time (yr)'), ylabel('v_1 (m)')
set(gca,'Fontsize',15)

tvals = (tresamp-tvec_eq(neq))./86400;% in days
T = table(tvals,u1);
% writetable(T,'assignment_week7/gps_timeseries.dat');
writetable(T,'gps_timeseries.dat');

















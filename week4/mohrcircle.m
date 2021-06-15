% script to Mohr's Circle for a given 2-d stress state
% Rishav Mallick, 2021, EOS

clc;
clear;

% inputs - provide an input stress field
% can be deviatoric or not, doesn't matter - we will plot it anyway

sx = -5;%input('sigma_x=');
sy = -1;%input('sigma_y=');
txy = -1;%input('tau_xy=');
theta = 0;%input('theta=');


%% if instead of stress tensor, this was a homogenous strain tensor, we could plot the dispalcements resulting from it
figure(10),clf
xval = linspace(-1,1,10)';
yval = linspace(-1,1,10)';
[xg,yg] = meshgrid(xval,yval);
ux = sx*xg + txy*yg;
uy = txy*xg + sy*yg;

plot(0,0,'ro','MarkerFaceColor','r'), hold on
quiver(xg,yg,ux,uy,'k','LineWidth',2)
axis tight equal
xlabel('x'),ylabel('y')
set(gca,'FontSize',15)

%% stress rotation by theta (if you want to)
disp('Stress Transformation')
sxp=((sx+sy)/2)+(((sx-sy)/2)*cosd(2*theta))+(txy*sind(2*theta));
syp=((sx+sy)/2)-((sx-sy)/2)*cosd(2*theta)-txy*sind(2*theta);
txyp=-((sx-sy)/2)*sind(2*theta)+txy*cosd(2*theta);
formatspec='sigma_x = %4.2f\nSigma_y= %4.2f\nTau_xy = %4.2f\n';
fprintf(formatspec,sxp,syp,txyp)

%% principle stresses
sig_avg=((sx+sy)/2);
R=sqrt(((sx-sy)/2)^2+txy^2);% max shear stress

sig1=sig_avg+R;
sig2=sig_avg-R;

% angles of principal stress and max shear
tp=.5*atan(2*txy/((sx-sy)))*180/pi;
ts=.5*atan(-(sx-sy)/(2*txy))*180/pi;
tp2=tp+90;
ts2=ts+90;

disp('Sigma Average and inplane Maximum Shear Stress');
formatspe='\nsigma_avg = %4.2f\nTau_max = %4.2f\n';
fprintf(formatspe,sig_avg,R);
disp('Angles Corresponding to Sigma_1');

xo=sig_avg;
formatsp='\nsigma_1 = %4.2f\nTheta_p1 = %4.2f\nTheta_s1 = %4.2f\n';
fprintf(formatsp,sig1,tp,ts);
formats='\nsigma_2 = %4.2f\nTheta_p2 = %4.2f\nTheta_s2 = %4.2f\n';
disp('Angles Sorresponding to Sigma_2');
fprintf(formats,sig2,tp2,ts2);

%% Mohr's circle
yo=0;
n=100;
r=R*ones(1,n);
theta1=linspace(0,2*pi,n);
[X,Y]=pol2cart(theta1,r);

% circle origin
X=X+xo;
Y=Y+yo;

x1=min([sx,sy]):1:max([sx,sy]);
y1=(txy/((((sx-sy)/2)+sig_avg)-sig_avg))*(x1-sig_avg);

figure(1),clf
plot(X,Y,'-','Color',rgb('steelblue'),'LineWidth',2), hold on
plot(x1,y1,'r-','LineWidth',1),
plot(0,0,'ks','MarkerFaceColor','r','MarkerSize',10)
set(gca,'Ydir','normal')
title('Mohrs Circle')
xlabel('\sigma_n \rightarrow')
ylabel('\tau \rightarrow')
hold on
grid on
axis equal
xC=sig_avg;
yC=0;
txtC=[' \sigma_{mean} = ',num2str(round(xC,1))];
text(xC,yC,txtC,'FontSize',15);

% label principal stresses
xb=sig2;
yb=0;
txt2= ['\leftarrow\sigma_{comp}= ',num2str(round(xb,1))];
text(xb,yb,txt2,'FontSize',15);
xa=sig1;
ya=0;
txt1=['\leftarrow\sigma_{ext}= ',num2str(round(xa,1))];
text(xa,ya,txt1,'FontSize',15);

% label the input stress state
txt3=['  \tau_{max}= ',num2str(round(R,1))];
plot(xC,R,'rp','MarkerFaceColor','r','MarkerSize',15)
plot([xC xC],[yC R],'ko--','LineWidth',1)
text(xC,R,{'',txt3},'FontSize',15);

x4=xC+((sx-sy)/2);
y4=txy;
txt4=['\sigma = ',num2str(x4)];
txt41=['\tau = ',num2str(y4)];
text(x4,y4,{txt4,txt41},'color','red','Fontsize',15);
x5=xC-((sx-sy)/2);
y5=(-txy);
txt5=['  \sigma = ',num2str(x5)];
txt51=['  \tau = ',num2str(y5)];
text(x5,y5,{txt5,txt51},'color','red','FontSize',15);
axis tight
set(gca,'FontSize',20)


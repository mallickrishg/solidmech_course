clear

ng = 40;
x = linspace(-5,5,ng)';
y = linspace(-5,5,ng)';
[X,Y] = meshgrid(x,y);

% temperature function
a = 2.5/pi;
T = sin(a*Y) + log(1+X.^2);

figure(1),clf
pcolor(x,y,T), hold on
axis tight equal, box on, grid on
xlabel('X'), ylabel('Y'), zlabel('T(X,Y)')
colormap(gray)

% compute gradients and plot it
[Tx,Ty] = gradient(T,x,y);
quiver(X,Y,Tx,Ty,'LineWidth',1)

% boat at (1,1,0) moving in the direction <1,1,1>
Tvalb = sin(a*1) + log(1+1^2);
% plot3(1,1,Tvalb,'rd','MarkerFaceColor','r','MarkerSize',10)
% quiver3(1,1,Tvalb,1,1,0,'r','LineWidth',2)
plot(1,1,'rd','MarkerFaceColor','r','MarkerSize',10)
quiver(1,1,3/2,1/2,'r','LineWidth',2)
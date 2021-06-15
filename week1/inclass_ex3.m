% calculating derivatives of a vector field
% Rishav Mallick, 2021, EOS

clear

% potential functions
u2 = @(m,a,x2,x3) m.*x2./(x2.^2 + x3.^2) + 2*a*x3.*cos(a.*x3.*x2);
u3 = @(m,a,x2,x3) m.*x3./(x2.^2 + x3.^2) + 2*a*x2.*cos(a.*x3.*x2);

m = 3;
a = 0.7;

% create x2,x3 vectors and a grid
ng = 200;

x2 = linspace(-3,3,ng)';
x3 = linspace(-3,3,ng);

[X2,X3] = meshgrid(x2,x3);

% compute u2 and u3 values 
u2_val = u2(m,a,X2,X3);
u3_val = u3(m,a,X2,X3);

% nan values outside the circle of r=3
rcond = sqrt(X2.^2 + X3.^2) > 3;
u2_val(rcond) = nan;
u3_val(rcond) = nan;

%% plot displacements
np = 2;

figure(1),clf
% plot the vector magnitude in color
contourf(x2,x3,sqrt(u2_val.^2 + u3_val.^2),[-5:5],'-','LineWidth',0.1), hold on
% plot vectors as arrows
quiver(X2(1:np:end),X3(1:np:end),u2_val(1:np:end),u3_val(1:np:end),'r','LineWidth',1)

axis tight equal
caxis([-1 1].*5)
cb=colorbar;cb.Label.String = '|u|';
xlabel('x_2'),ylabel('x_3')
set(gca,'FontSize',20,'LineWidth',1)
colormap(gray)

%% compute gradients for u2 and u3 separately
% fill this section out by yourself
% [u22,u23] = gradient(u2.....)
% [u32,u33] = gradient(u3.....)
% Plot the gradients of u2 and u3 separately
% Compute and plot the scalar field u22 + u33. What does it look like
% compared to the vector field (u2,u3)?
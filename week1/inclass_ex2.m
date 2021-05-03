% Script for in-class exercise 
% Rishav Mallick, EOS, 2021
clear

t = linspace(1e-2,2,1000)';

f = (6*(t.^(1/3)))./(4*t + 1);
fd = (-16*(t.^(1/3)) + 2./(t.^(2/3)))./((4*t+1).^2);
fd_n = gradient(f)./gradient(t);

figure(1),clf
plot(t,f,'LineWidth',2)
axis tight, grid on

figure(2),clf
plot(t,fd,'r-','LineWidth',2), hold on
plot(t,fd_n,'b--','LineWidth',2)
plot(t,0.*t,'k-','LineWidth',1)
axis tight, grid on

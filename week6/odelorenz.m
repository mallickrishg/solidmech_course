function Yp = odelorenz(~,Y,evl)

x = Y(1);
y = Y(2);
z = Y(3);

dxdt = evl.sigma*(y-x);
dydt = evl.r*x - z*x - y;
dzdt = x*y - evl.b*z;

Yp = zeros(size(Y));

Yp(1) = dxdt;
Yp(2) = dydt;
Yp(3) = dzdt;

end
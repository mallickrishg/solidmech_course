function Yp = ode_rsf(~,Y,evl)

damping = 0.5*evl.K;

% State variables
th = Y(2);

% Slip velocities 
V = evl.Vo.*exp(Y(3));

% Initiate state derivative
Yp=zeros(size(Y));  

% Slip velocity
Yp(1) = V;

% Rate of state (rate of log(theta/theta0))
dth=(evl.Vo.*exp(-th)-V)./evl.l;

Yp(2) = dth;

kv = -evl.K*(V - evl.Vpl);

% stress rate
Yp(3) = (kv - evl.b.*evl.sigma.*dth)./(evl.a.*evl.sigma + damping.*V);

end
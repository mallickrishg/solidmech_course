function Yp = odeElastic_a(~,Y,ss,evl)
% ODE function for faults and shear zones evolution
% Y = [slip,theta,V... - flt
% We express dY/dt as f(Y) and using RK-4th order to integrate with
% adaptive timesteps
% Rishav Mallick, 2019, EOS

%% FAULTS
G = 30e3; % shear modulus
damping = 0.5*G./ss.Vs;

% State variables
th   = Y(2:ss.dgf:ss.M*ss.dgf);

% Slip velocities 
V    = ss.Vo.*exp(Y(3:ss.dgf:ss.M*ss.dgf));

% Initiate state derivative
Yp=zeros(size(Y));  

% Slip velocity
Yp(1:ss.dgf:ss.M*ss.dgf) = V;

% Rate of state (rate of log(theta/theta0))
dth=(ss.Vo.*exp(-th)-V)./ss.l;
Yp(2:ss.dgf:ss.M*ss.dgf) = dth;

%% INTERACTIONS
% Acceleration (rate of log(V/Vo))
kv = evl.K*V + evl.tau0;
Yp(3:ss.dgf:ss.M*ss.dgf) = (kv - ss.b.*ss.sigma.*dth)./(ss.a.*ss.sigma + damping.*V);

end
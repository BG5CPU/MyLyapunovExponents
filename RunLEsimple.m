clear;clc;
% x(1) = id
% x(2) = iq
% x(3) = omegaP
xInit = [1,1,1];
orthmatrix = [1 0 0;
            0 1 0;
            0 0 1];
% GAMA = parameter(1); 
% SIGMA = parameter(2);
parameter = [0.15, 0.20, 10.0];
% u(1) = ud
% u(2) = uq
% u(3) = TL
u = [0, 0, 0];
x0 = zeros(12,1);
x0(1:3) = xInit;
x0(4:12) = orthmatrix;
tstart = 0; % start time
tstep = 1e-3; % step time
wholetimes = 1e5; % total number of cycles
steps = 10; % steps per evolution
iteratetimes = wholetimes/steps; % numbers of evolution
mod = zeros(3,1);
lp = zeros(3,1);
Lyapunov1 = zeros(iteratetimes,1);
Lyapunov2 = zeros(iteratetimes,1);
Lyapunov3 = zeros(iteratetimes,1);
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);

for i=1:iteratetimes
    tspan = tstart:tstep:(tstart + tstep*steps);    
    [t, x] = ode45(@(t, x)Rossler_ly(t, x, u, parameter), tspan, x0, opts);
    % value of the last instant of integration
    x0 = x(size(x,1),:);
    % re-define the start time
    tstart = tstart + tstep*steps;
    xx0 = [x0(4) x0(7) x0(10);
        x0(5) x0(8) x0(11);
        x0(6) x0(9) x0(12)];
    % Orthogonalization
    xx0 = ThreeGS(xx0);
    % vector module
    mod(1) = sqrt(xx0(:,1)'*xx0(:,1));
    mod(2) = sqrt(xx0(:,2)'*xx0(:,2));
    mod(3) = sqrt(xx0(:,3)'*xx0(:,3));
    xx0(:,1) = xx0(:,1)/mod(1);
    xx0(:,2) = xx0(:,2)/mod(2);
    xx0(:,3) = xx0(:,3)/mod(3);
    lp = lp+log(abs(mod));
    %三个Lyapunov指数
    Lyapunov1(i) = lp(1)/(tstart);
    Lyapunov2(i) = lp(2)/(tstart);
    Lyapunov3(i) = lp(3)/(tstart);
    x0(4:12) = xx0';  
end

i = 1:iteratetimes;
plot(i,Lyapunov1,i,Lyapunov2,i,Lyapunov3)

















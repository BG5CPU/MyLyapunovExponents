function [OUT] = PmsmLorenzfuncU(t, x, parameter, tud, ud)
% x(1) = id
% x(2) = iq
% x(3) = omegaP
% x(4)~~x(12) Jcobian
% u(1) = ud
% u(2) = uq
% u(3) = TL

GAMA = parameter(1); 
SIGMA = parameter(2);

ud = interp1(tud, ud, t);

Q = [x(4), x(7), x(10);
    x(5), x(8), x(11);
    x(6), x(9), x(12)];

%Lorenz equation
dx1 = -x(1)+x(2)*x(3)+ud;
dx2 = -x(2)-x(1)*x(3)+GAMA*x(3)+0;
dx3 = SIGMA*(x(2)-x(3))+0;

%Linearized system
JcoM = [-1, x(3), x(2);
        -x(3), -1, -x(1)+GAMA;
        0, SIGMA, -SIGMA];

%Variational equation
F = JcoM*Q;

%Output data must be a column vector
OUT=[dx1; dx2; dx3; F(:)];

end


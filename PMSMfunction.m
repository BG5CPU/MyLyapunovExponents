function [dx] = PMSMfunction(t, x, u, parameter)
% x(1) = id
% x(2) = iq
% x(3) = omegaP
% u(1) = ud
% u(2) = uq
% u(3) = TL
% parameter(1) = Rs % 0.011
% parameter(2) = Ld % 0.018
% parameter(3) = Lq % 0.018
% parameter(4) = Fi % 0.055
% parameter(5) = p % 4
% parameter(6) = J % 0.1
% parameter(7) = b % 0.2
dx = [
    -(parameter(1)/parameter(2))*x(1)+(parameter(3)/parameter(2))*x(2)*x(3)+(1/parameter(2))*u(1);
%    -(parameter(1)/parameter(3))*x(2)-(parameter(2)/parameter(3))*x(1)*x(3)-(parameter(4)/parameter(3))*x(3)+(1/parameter(3))*u(2);
    -(parameter(1)/parameter(3))*x(2)-(parameter(2)/parameter(3))*x(1)*x(3)+(parameter(4)/parameter(3))*x(3)+(1/parameter(3))*u(2);
    (parameter(5)^2*parameter(4)/parameter(6))*x(2)-(parameter(7)/parameter(6))*x(3)-(parameter(5)/parameter(6))*u(3);
];
end


function [dxdt] = PMSMfunFlex(t, x, parameter)

% parameter = [A1,A2,A3,A4,C1,C2,C3,C4,C5,D1,D2,D3,D4,D5,D6,D7,D8];
A1 = parameter(1);
A2 = parameter(2);
A3 = parameter(3);
A4 = parameter(4);
C1 = parameter(5);
C2 = parameter(6);
C3 = parameter(7);
C4 = parameter(8);
C5 = parameter(9);
D1 = parameter(10);
D2 = parameter(11);
D3 = parameter(12);
D4 = parameter(13);
D5 = parameter(14);
D6 = parameter(15);
D7 = parameter(16);
D8 = parameter(17);

% ud = interp1(tud, ud, t);
ud = 0.1375*sin(2*pi*t);

Q = [x(6), x(11), x(16), x(21), x(26);
     x(7), x(12), x(17), x(22), x(27);
     x(8), x(13), x(18), x(23), x(28);
     x(9), x(14), x(19), x(24), x(29);
     x(10),x(15), x(20), x(25), x(30)];

% system function
dx1 = x(2);
dx2 = x(3);
dx3 = A1*x(1)+A2*x(2)+A3*x(3)+A4*x(5);
dx4 = C1*x(4)+C2*x(1)*x(5)+C3*x(2)*x(5)+C4*x(3)*x(5)+C5*ud;
dx5 = D1*x(5)+D2*x(1)*x(4)+D3*x(2)*x(4)+D4*x(3)*x(4)+D5*x(1)+D6*x(2)+D7*x(3)+D8*0;

% Linearized system
JcoM = [0, 1, 0, 0, 0;
        0, 0, 1, 0, 0;
        A1, A2, A3, 0, A4
        C2*x(5), C3*x(5), C4*x(5), C1, C2*x(1)+C3*x(2)+C4*x(3);
        D2*x(4)+D5, D3*x(4)+D6, D4*x(4)+D7, D2*x(1)+D3*x(2)+D4*x(3), D1];

% Variational equation
F = JcoM*Q;

% Output data must be a column vector
dxdt = [dx1; dx2; dx3; dx4; dx5; F(:)];

end












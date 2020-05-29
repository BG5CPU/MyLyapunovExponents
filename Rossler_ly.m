function dX = Rossler_ly(t,X,u,parameter)
%  Rossler�����ӣ���������Lyapunovָ��
%       a=0.15,b=0.20,c=10.0
%       dx/dt = -y-z,
%       dy/dt = x+ay,
%       dz/dt = b+z(x-c),
% a = 0.15;
% b = 0.20;
% c = 10.0;
a = parameter(1);
b = parameter(2);
c = parameter(3);
x=X(1); y=X(2); z=X(3);
% Y������������Ϊ�໥�����ĵ�λ����
Y = [X(4), X(7), X(10);
    X(5), X(8), X(11);
    X(6), X(9), X(12)];
% ��������ĳ�ʼ�����ز�����
dX = zeros(12,1);
% Rossler������
dX(1) = -y-z+u(1);
dX(2) = x+a*y+u(2);
dX(3) = b+z*(x-c)+u(3);
% Rossler�����ӵ�Jacobi����
Jaco = [0 -1 -1;
       1 a  0;
       z 0  x-c];
dX(4:12) = Jaco*Y;








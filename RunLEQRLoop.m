clear;clc;
% ItrNum=10;
linODEnum = 9; %No. of linearized ODEs
dim = sqrt(linODEnum); %Dimension of the linearized system
% x(1) = id
% x(2) = iq
% x(3) = omegaP
xInit = [1,1,1];
Q0 = eye(dim);
% GAMA = parameter(1); 
% SIGMA = parameter(2);
parameter = [24.4444 2.8877];
% u(1) = ud
% u(2) = uq
% u(3) = TL
u = [0, 0, 0];

%Initial conditions
IC0 = zeros(12,1);
IC0(1:3) = xInit;
IC0(4:12) = Q0;

DiscardItr = 200; %Iterations to be discarded
                  % Transient iterations to be discarded:
                  % 200 Iterations = 200*10 time steps = 20 sec
InitialTime = 0; %Initial time
FinalTime = 400; %Final Time: 1000 sec
TimeStep = 1e-2; % time step
UpdateStepNum = 10; %Lyapunov exponents updating steps
Iteration = UpdateStepNum*TimeStep; % 10*0.01=0.1
DiscardTime = DiscardItr*Iteration+InitialTime; % 200*0.1=20

tud = linspace(0,FinalTime,FinalTime/TimeStep);
ud = 100*sin(2*pi*tud);

opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
LyValue = [];

for GAMA = 1:1:50
    parameter(1) = GAMA;

    T1 = InitialTime;
    T2 = T1+Iteration;
    TSpan = [T1:TimeStep:T2];
    
    n = 0; % Iteration counter
    k = 0; % Effective iteration counter
         % (discarded iterations are not counted)
    h = 1;
    Sum = zeros(1,dim);
    xData = [];
    yData = [];
    A = [];
    IC = IC0;
    
    %Main loop
    while (1)
        n = n+1;
        %Integration
        [t, x] = ode45(@(t, x)PmsmLorenzfuncU(t, x, parameter, tud, ud), TSpan, IC, opts);
        [rX,cX]=size(x);
        L=cX-linODEnum; %No. of initial conditions for
                        %the original system
        for i=1:dim
            m1=L+1+(i-1)*dim;
            m2=m1+dim-1;
            A(:,i)=(x(rX,m1:m2))';
        end
        %QR decomposition
        [Q,R]=qr(A);
        if T2>DiscardTime
            Q0=Q;
         else
            Q0=eye(dim);
        end
        %Any zero diagonal element will cause overflow
        %in the following calculation, so discard this step.
        permission=1;
        for i=1:dim
            if R(i,i) == 0
                permission = 0;
                break;
            end
        end
        %To determine the Lyapunov exponents
        if (T2>DiscardTime & permission)
            k=k+1;
            T=k*Iteration;
            TT=n*Iteration+InitialTime;
            %There are d Lyapunov exponents
            Sum=Sum+log(abs(diag(R))');
            lambda=Sum/T;
            %Sort the Lyapunov exponents in descenting order
            Lambda=fliplr(sort(lambda));
            %To calculate the Lyapunov dimension (or Kaplan-Yorke dimension)
            LESum=Lambda(1);
            LD=0;
            if Lambda(1)>0
                for N=1:dim-1
                    if Lambda(N+1)~=0
                        LD=N+LESum/abs(Lambda(N+1));
                        LESum=LESum+Lambda(N+1);
                        if LESum<0
                            break;
                        end
                    end
                end
            end        
            xData=[xData;TT];
            yData=[yData;lambda];        
        end
    
        %If calculation is finished, exit the loop.
        if T2>=FinalTime
            break;
        end
    
        %Update the initial conditions and time span for the next iteration
        xInit = x(rX,1:L);
        T1=T1+Iteration;
        T2=T2+Iteration;
        TSpan=[T1:TimeStep:T2];
        IC=[xInit(:);Q0(:)];
    end

    LyValue(1,length(LyValue)+1) = GAMA;
    LyValue(2:4,length(LyValue)) = yData(size(yData,1),:);

end

Zeroo = zeros(1,length(LyValue));

figure('color',[1 1 1]);
set(gcf,'position',[50 50 500 400]);
plot(LyValue(1,:), Zeroo, ':', 'LineWidth', 1, 'color', [0, 0, 0]);
hold on
plot(LyValue(1,:), LyValue(2,:), '-', 'LineWidth', 2, 'color', [0, 0, 1]);
hold on
plot(LyValue(1,:), LyValue(3,:), '-', 'LineWidth', 1.5, 'color', [1, 0, 0]);
hold on
plot(LyValue(1,:), LyValue(4,:), '-', 'LineWidth', 1.5, 'color', [0, 1, 0]);
set(gca,'XLim',[1 25]);
legend('0 line', '\lambda_1','\lambda_2','\lambda_3');
set(gca,'fontsize',15,'fontname','Times New Roman');
xlabel('\gamma','Fontname', 'Times New Roman','FontSize',16,'FontAngle','italic');
ylabel('Lyapunov Exponents','Fontname', 'Times New Roman','FontSize',16);
grid on;
set(gca,'position',[0.13 0.15 0.83 0.82]);

















clear;clc;
% parameter
Rs = 0.011;
Ld = 0.018;
Fi = 0.055;
p = 4;
JM = 0.002; JL = 0.1;
BM = 0.18; BS = 0.002;
KS = 500;
% ItrNum=10;
linODEnum = 25; %No. of linearized ODEs
dim = sqrt(linODEnum); %Dimension of the linearized system
% x(1) = id
% x(2) = iq
% x(3) = omegaP
xInit = [0,0,0,1,1];
Q0 = eye(dim);
% u(1) = ud
% u(2) = uq
u = [0, 0];

%Initial conditions
IC0 = zeros(30,1);
IC0(1:5) = xInit;
IC0(6:30) = Q0;

DiscardItr = 200; %Iterations to be discarded
                  % Transient iterations to be discarded:
                  % 200 Iterations = 200*10 time steps = 20 sec
InitialTime = 0; %Initial time
FinalTime = 400; %Final Time: 1000 sec
TimeStep = 1e-2; % time step
UpdateStepNum = 10; %Lyapunov exponents updating steps
Iteration=UpdateStepNum*TimeStep; % 10*0.01=0.1
DiscardTime=DiscardItr*Iteration+InitialTime; % 200*0.1=20

T1=InitialTime;
T2=T1+Iteration;
TSpan=[T1:TimeStep:T2];

opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
LyValue = [];

for para = 10:10:10000
    KS = para;
    
    T1 = InitialTime;
    T2 = T1+Iteration;
    TSpan = [T1:TimeStep:T2];
    
    n=0; % Iteration counter
    k=0; % Effective iteration counter
         % (discarded iterations are not counted)
    h=1;
    Sum=zeros(1,dim);
    xData=[];
    yData=[];
    A=[];
    IC = IC0;

    % calculate the paremeter
    A1 = -(BM*KS)/(JM*JL);
    A2 = -(JM*KS+JL*KS+BM*BS)/(JM*JL);
    A3 = -(JL*BS+JM*BS+JL*BM)/(JM*JL);
    A4 = (p*Fi)/(JM*JL);
    C1 = -Rs/Ld; C2 = KS*p; C3 = BS*p; C4 = JL*p; C5 = 1/Ld;
    D1 = -Rs/Ld; D2 = -KS*p; D3 = -BS*p; D4 = -JL*p;
    D5 = (Fi*KS*p)/Ld; D6 = (Fi*BS*p)/Ld;
    D7 = (Fi*JL*p)/Ld; D8 = 1/Ld;
    parameter = [A1,A2,A3,A4,C1,C2,C3,C4,C5,D1,D2,D3,D4,D5,D6,D7,D8];
    B3 = JL; B2 = BS; B1 = KS;

    %Main loop
    while (1)
        n = n+1;
        %Integration
        [t, x] = ode45(@(t, x)PMSMfunFlex(t, x, parameter), TSpan, IC, opts);
        [rX,cX]=size(x);
        LL=cX-linODEnum; %No. of initial conditions for
                        %the original system
        for i=1:dim
            m1=LL+1+(i-1)*dim;
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
        xInit = x(rX,1:LL);
        T1=T1+Iteration;
        T2=T2+Iteration;
        TSpan=[T1:TimeStep:T2];
        IC=[xInit(:);Q0(:)];
    end
    LyValue(1,length(LyValue)+1) = para;
    LyValue(2:6,length(LyValue)) = yData(size(yData,1),:);
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
hold on
plot(LyValue(1,:), LyValue(5,:), ':', 'LineWidth', 4, 'color', [0.06275, 0.30588, 0.5451]);
hold on
plot(LyValue(1,:), LyValue(6,:), '-', 'LineWidth', 1.5, 'color', [0.7451,0.7451,0.7451]);
set(gca,'XLim',[0.001 0.5]);
legend('0 line', '\lambda_1','\lambda_2','\lambda_3','\lambda_4','\lambda_5');
set(gca,'fontsize',15,'fontname','Times New Roman');
xlabel('$B_S$','interpreter','latex','Fontname', 'Times New Roman','FontSize',16,'FontAngle','italic');
ylabel('Lyapunov Exponents','Fontname', 'Times New Roman','FontSize',16);
grid on;
set(gca,'position',[0.14 0.15 0.825 0.82]);


Zeroo = zeros(1,length(LyValue));
figure('color',[1 1 1]);
set(gcf,'position',[50 50 500 400]);
plot(LyValue(1,:), Zeroo, ':', 'LineWidth', 1, 'color', [0, 0, 0]);
hold on
plot(LyValue(1,:), LyValue(2,:), '-', 'LineWidth', 3, 'color', [0, 0, 1]);
hold on
plot(LyValue(1,:), LyValue(3,:), '-', 'LineWidth', 3, 'color', [1, 0, 0]);
hold on
plot(LyValue(1,:), LyValue(4,:), '-', 'LineWidth', 1.5, 'color', [0, 1, 0]);
hold on
plot(LyValue(1,:), LyValue(5,:), ':', 'LineWidth', 4, 'color', [0.06275, 0.30588, 0.5451]);
hold on
plot(LyValue(1,:), LyValue(6,:), '-', 'LineWidth', 1.5, 'color', [0.7451,0.7451,0.7451]);
set(gca,'XLim',[0.1 0.3],'YLim',[-2 1.5]);
set(gca,'fontsize',15,'fontname','Times New Roman');
grid on;
set(gca,'position',[0.12 0.15 0.825 0.82]);

















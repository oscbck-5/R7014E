h1 = 2;
h2 = 2;
h3 = 2;
h4 =2;

A1 = 28; A3 = A1;
A2 = 32; A4 = A2;
a1 = 0.071; a3 = a1;
a2 = 0.057; a4 = a2;
g  = 981;
kc = 0.5;
k1 = 3.33;
k2 = 3.35;
gamma1 = 0.7;
gamma2= 0.6;




T1 = A1/a1*sqrt(2*h1/g);
T2 = A2/a2*sqrt(2*h2/g);
T3 = A3/a3*sqrt(2*h3/g);
T4 = A4/a4*sqrt(2*h4/g);

c1 =T1*k1*kc/A1; 
c2= T2*k2*kc/A2;
lambda1 = 10;
lambda2 = 5;


G1 =[gamma1*c1/(T1*s+1)];
G2 = [gamma2*c2/(T2*s+1)];
G12 = [(1-gamma2)*c1/((T1*s+1)*(T3*s+1))] ;
G21 = [(1-gamma1)*c2/((T2*s+1)*(T4*s+1))] ;
G = [G1 G12; G21 G2];
%Fy11= eye(1)*(T1*s+1)/((gamma1*c1)*(lambda*s+1)) - eye(1)/(lambda*s+1)*(T1*s+1)/(gamma1*c1*(lambda*s+1))

Q1 = G1\(1/(lambda1*s+1));
Fy1 = (eye(1)-Q1*G1)\Q1;

Q2 = G2\(1/(lambda2*s+1));
Fy2 = (eye(1)-Q2*G2)\Q2;

Fy1 = minreal(Fy1)
Fy2 = minreal(Fy2)
%Q2 = 1/(lambda*s+1)

[num den] = tfdata(Fy1(1,1), 'v')
[num2 den2] = tfdata(Fy2(1,1), 'v')



%% Kalman

states = minreal(ss(G))
[A,B,C,D] =ssdata(states)

%Kalman.
%R= [1 0;  0 1];
Qk = 1*eye(4);

[P,L,K] =care(A',C',Qk)   %ehm?
%K = K;

Ltest1 = eye(2)*B'*P
Lrtest1 = (C*(B*Ltest1-A)^-1*B)^-1
%K = (P*C'*R(2,2))


%% L statefeedback


%% integral ACTION
clc
Ai = [A zeros(4,2);-C zeros(2,2)]
Bi = [B zeros(4,2); zeros(2,2) eye(2)]
Ci = [C zeros(2,2); zeros(2,4) eye(2)]

%M = Ci;
M = [-C zeros(2,2)];
%Lr = Bi*Ltest'

%Q1 = eye(6);
%Q2 = eye(2);
% Q1 = [1 0   0   0   0   0;
%          0  1   0   0   0   0;
%          0  0   1   0   0   0;
%          0  0   0   1   0   0;
%          0  0   0   0   1   0;
%          0  0   0   0   0   1];
Q1 = diag([1 0.05 1 1 1 0.01]); %FAST h1, neglect h2
Q2 = diag([1 0.01]);                  %FAST h1, neglect h2

Q1 = diag([5 8 1 1 5 8]);          %Both equally fast
Q2 = diag([20 80]);                  %Both equally fast


[P2,L2,K2] =care(Ai,Bi(1:end,1:2),Q1,Q2)



L = eye(2)*Bi(1:end,1:2)'*P2
sim('R7014E_QuadTankSim_With_Kalman.slx',1000)
%Lr = (M*(Bi*Ltest-Ai)^-1*Bi(1:end,3:4 ))^-1
%Lr = (M*(Bi(1:end,1:2)'*Ltest-Ai)^-1*Bi(1:end,1:2))

%L = Ltest(1:2,1:4)  %kaos


%Lr = (M*((Bi*Ltest-A)^-1)*Bi)^-1;

%% 
clc
[num2,den2] = tfdata(G) 
celldisp(num2)
celldisp(den2)











%%
clc

%M = Ci;
M = [C zeros(2,2)];
%Lr = Bi*Ltest'

Q2 = eye(2);
[P2,L2,K2] =care(Ai',Ci',Q1(1,1)*eye(6),Q1)
%Ltest = eye(4)*Bi'*P2

%Lr = (M*(Bi*Ltest-Ai)^-1*Bi(1:end,3:4 ))^-1


%L = Ltest(1:2,1:4)  %kaos















%%

Q2 = eye(6);
clc
[P2,L2,K2] =care(Ai,Bi,Q2)
M = C;

%Lr = (M*((Bi*Ltest-A)^-1)*Bi)^-1;
Lrtest = ((M*(B*Ltest-A)^-1)*B)^-1


%[P2,L2,K2] =care(A,B,eye(4)*Q1(1,1))
%Ltest = 1.*inv(eye(2))*B'*P2
%M = C; 


%Lr = (M*((B*Ltest-A)^-1)*B)^-1;

%K=(P*C2'+R(1,2))*inv(R(2,2))

%%
%Rmse
% imse(height.signals(1).values(120:end))
height1_mean = mean((10-height.signals(1).values(120:end)).^2)
RMSE1 = sqrt(height1_mean)

height2_mean = mean((10-height.signals(2).values(120:end)).^2)
RMSE2 = sqrt(height2_mean)
%% 












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
gamma2 = 0.6;




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
%%

states = minreal(ss(G))
[A,B,C,D] =ssdata(states)

%Kalman.
R= [1 0
   0 1];

[X,L,K] =care(A',C',eye(4)*R(1,1))   %ehm?
K = K';

Lr = (M*((B*K'-A)^-1)*B)^-1

%K=(P*C2'+R(1,2))*inv(R(2,2))

M = C; 

















%% System Param
%h vs height Voltage in cm
clc
h = 5:1:15;
hV = [1.826, 2.539, 3.452, 4.18, 4.99, 5.762, 6.665, 7.36, 8.193, 8.93, 9.66];

%stem(h,hV) %linear

%Voltage vs time for 0-15cm
v = 1:1:10;
vt = [21.18, 26.66 30.97, 34.28, 37.92, 39.73, 42.53, 45.45, 46.49, 47.67];
%plot(v,vt) % linear njaaaa

%Kalman filter
C1=[1 0];
C2=[0 1];
C=[C1;C2];
R= [0.1 0
   0 0.1];

g = 981;
%A = 3.5e-3;
A=33;
%a = 0.16e-4;
a=0.17;
h10=7;
h20=h10;
h0=0;
Ac=[-(a*g*sqrt(2))/(2*sqrt((h10+h0)*g)*A) 0
    (a*g*sqrt(2))/(2*sqrt((h10+h0)*g)*A) -(a*g*sqrt(2))/(2*sqrt((h20+h0)*g)*A)];

Bc=[1/A
    0];
N=[1
   0];

[P L K]=care(Ac',C2',eye(2)*R(1,1),R(2,2))
K = K';
%K=(P*C2'+R(1,2))*inv(R(2,2))    
%% Discrete
clc
Plant=ss(Ac,Bc,C2,0);
%[kalmf,K,P,M] = kalman(Plant,1,1);


PlantD=c2d(Plant,0.005);
[Ak,Bk,Ck,Dk]=ssdata(PlantD);

[P L K]=dare(Ak',Ck',eye(2)*R(1,1),R(2,2))
K = K'
%K = K'
%K=(P*Ck'+R(1,2))*inv(R(2,2))



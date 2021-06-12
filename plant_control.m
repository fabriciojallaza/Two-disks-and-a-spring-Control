%% Plant definition
%Variables
th1=0.0019;
th2=0.0019;
kt1=3.6242*exp(-005);
kt2=2.05*exp(-005);
c=0.2038;
d=0.0022;
d1=4.7*exp(005);
d2=4.7*exp(005);

%Plant
A=[0 1 0 0;-c/th1 -(d1+d)/th1 -c/th1 -d/th1;0 0 0 1; -c/th2 -d/th2 -c/th2 -(d2+d)/th2];
B=[0 0;kt1/th1 0; 0 0; 0 kt2/th2];
C=[1 0 0 0; 0 0 1 0];
D=[0 0;0 0];
plant=ss(A,B,C,D);

%% Controllability and observability
% MC=[B A*B A^2*B];
% MO=[C;C*A;C*A^2];
% cont=det(MC);
% obs=det(MO);

%% Requirements
PO=6; %overshoot less than 6%
rt=2; %rise tiem less than 2 seconds
st=5; %settling time less than 5 seconds
sse=2; %Steady-state error less than 2%

%% Control implementation
%Determination of the poles
p=pi;
e=(abs(log(PO/100)))/(sqrt(p^2+(log(PO/100))^2));
wn_prima=4/(e*st);
p1=-e*wn_prima+wn_prima*sqrt(e^2-1);
p2=-e*wn_prima-wn_prima*sqrt(e^2-1);
p3=real(p1)*10;
p4=real(p1)*10-1;

%Control implementation
k=place(A,B,[p1 p2 p3,p4])
plant_uc=ss(A-B*k, B, C-D*k, D);

%% LQR
Q =  [19000   0     0     0;
      0     0     0     0;
      0     0     19000     0;
      0     0     0     0];
 
R = 0.00001;

[k1,S,e] = lqr(A,B,Q,R);
Kpre=inv(-(C-D*k1)*inv(A-B*k1)*B+D);

plant_uc2=ss(A-B*k1, B, C, D);
% 
[y1,t]=step(plant_uc2*Kpre); grid on
stepResults=stepinfo(y1,t);
aux2=stepResults.SettlingTime
aux3=stepResults.Overshoot

%% Simulations
subplot(311);
step(plant), grid on, title('PLANT')
subplot(312);
step(plant_uc), grid on, title('PLANT UNDER CONTROL')
subplot(313);
step(plant_uc2), grid on, title('PLANT UNDER CONTROL')


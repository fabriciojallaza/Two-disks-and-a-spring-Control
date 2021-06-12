%% Plant definition
%Variables
th1=0.0019;
th2=0.0019;
kt1=3.6242e-005;
kt2=2.05e-005
c=0.2038;
d=0.0022;
d1=4.7e-005;
d2=4.7e-005;

%Plant
A=[0 1 0 0;-c/th1 -(d1+d)/th1 -c/th1 -d/th1;0 0 0 1; -c/th2 -d/th2 -c/th2 -(d2+d)/th2];
B=[0 0;kt1/th1 0; 0 0; 0 kt2/th2];
C=[1 0 0 0; 0 0 1 0];
D=[0 0;0 0];
plant=ss(A,B,C,D);

%% Controllability and observability
Co = ctrb(A,B);
Ob = obsv(A,C);
rc = rank(Co)
ro = rank(Ob)

%% Requirements
PO=6; %overshoot less than 6%
rt=2; %rise tiem less than 2 seconds
st=5; %settling time less than 5 seconds
sse=2; %Steady-state error less than 2%

%% Control implementation
%Pole placement
p=pi;
E=(abs(log(PO/100)))/(sqrt(p^2+(log(PO/100))^2));
wn_prima=4/(E*st);
p1=-E*wn_prima+wn_prima*sqrt(E^2-1);
p2=-E*wn_prima-wn_prima*sqrt(E^2-1);
p3=real(p1)*10;
p4=real(p1)*100;

k1=place(A,B,[p1 p2 p3,p4])
plant_uc=ss(A-B*k1, B, C-D*k1, D);

%% Static error
inv(A-B*k1)*B-D
%% 
%Ganancia Kpre
kpre=1/-((C-D*k1)*inv(A-B*k1)*B-D);
plant_uckpre=ss(A-B*k1,B*kpre,C-D*k1,D*kpre);

%Con parte Integral
n=size(A);
n=n(1);
m=size(C);
m=m(1);
Aext=[A zeros(n,1); -C zeros(m,1)];
Bext=[B; -D];
p5=real(p1)*10-1;
kext=place(Aext,Bext,[p1,p2,p3,p4,p5]);
kn=kext(1:n);
ke=kext(n+1);
plant_ucint=ss([A-B*kn -B*ke;-C+D*kn D*ke],[zeros(n,1); 1],[C 0],0);

%% LQR
Q =  [19000   0     0     0;
      0     0     0     0;
      0     0     19000     0;
      0     0     0     0];
 
R = 0.00001;

[k1,S,E] = lqr(A,B,Q,R);
Kpre=inv(-(C-D*k1)*inv(A-B*k1)*B+D);

plant_uc2=ss(A-B*k1, B, C, D);
% 
[y1,t]=step(plant_uc2*Kpre); grid on
stepResults=stepinfo(y1,t);
aux2=stepResults.SettlingTime
aux3=stepResults.Overshoot

%% ITAE criteria
%Determination of the poles
itaeeq1=[1 2.1*wi 3.4*wi^2 2.7*wi^3 wi^4];
itaepol1=roots(itaeeq1);
pi1=itaepol1(1);
pi2=itaepol1(2);
pi3=itaepol1(3);
pi4=itaepol1(4);
itaeeq2=[1 2.8*wi 5.0*wi^2 5.5*wi^3 +3.4*wi^4 wi^5];
itaepol2=roots(itaeeq2);
pin1=itaepol2(1);
pin2=itaepol2(2);
pin3=itaepol2(3);
pin4=itaepol2(4);
pin5=itaepol2(5);

%Control implementation
ki=place(A,B,[pi1,pi2,pi3,pi4]);
PlantITAE=ss(A-B*ki, B, C-D*ki, D);

%%Static error
%Kpre
kprei=1/((C)*(inv(-A+B*ki))*B);
PlantITAEk=ss(A-B*ki,B*kprei,C-D*ki,D*kprei);
%Integral
% kiext=place(Aext,Bext,[pin1,pin2,pin3,pin4]);
% ki2=kiext(1:n);
% kei=kiext(n+1);
% PlantITAEi=ss([A-B*ki2 -B*kei;-C+D*ki2 D*kei],[zeros(n,1); 1],[C 0],0);

%% Simulations
subplot(311);
step(plant), grid on, title('PLANT')
subplot(312);
step(plant_uc), grid on, title('PLANT UNDER CONTROL')
subplot(313);
step(plant_uc2), grid on, title('PLANT UNDER CONTROL')


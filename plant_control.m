clc
clear
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
B=[0 ;kt1/th1 ; 0 ; 0 ]; % we are interested in the first motor torque
C=[ 0 0 1 0]; % we are intersted in the second motor output/angle
D=[0];
plant=ss(A,B,C,D);

%% Controllability and observability
Co = ctrb(A,B);
Ob = obsv(A,C);
Rango_Controlablidad = rank(Co)
Rango_Observabilidad = rank(Ob)
%With a rank of 4 each we see that are controlable and observable
%% Requirements
OS=6; %overshoot less than 6%
rt=2; %rise tiem less than 2 seconds
st=5; %settling time less than 5 seconds
sse=2; %Steady-state error less than 2%

%% Control implementation
%Pole placement
p=pi;
E=(abs(log(OS/100)))/(sqrt(p^2+(log(OS/100))^2));
wn_prima=4/(E*st);
p1=-E*wn_prima+wn_prima*sqrt(E^2-1);
p2=-E*wn_prima-wn_prima*sqrt(E^2-1);
p3=real(p1)*10;
p4=real(p1)*100;
p=[p1 p2 p3 p4];
k1=place(A,B,[p1 p2 p3,p4]);
plant_uc=ss(A-B*k1, B, C-D*k1, D); %with pre compensation gai

%Pole placement with statit error
%Kpre Gain
% kpre=1/-((C-D*k1)*inv(A-B*k1)*B-D);
% plant_uckpre=ss(A-B*k1,B*kpre,C-D*k1,D*kpre);

%Integral part
n=size(A);
n=n(1);
m=size(C);
m=m(1);
% a, b extended matrix
Aext=[A zeros(n,1); -C zeros(m,1)];
Bext=[B; -D];
% Contrabilidad Matrices extendidas
Coext = ctrb(Aext,Bext);
Rango_Controlablidad_Mext = rank(Coext)

p5=real(p1)*10-1;
kext=place(Aext,Bext,[p1,p2,p3,p4,p5]);
kn=kext(1:n);
ke=kext(n+1);
plant_uci=ss([A-B*kn -B*ke;-C+D*kn D*ke],[zeros(n,1); 1],[C 0],0);

%% ITAE criteria
%Determination of the poles
wi=1.5;
itaeeq1=[1 2.1*wi 3.4*wi^2 2.7*wi^3 wi^4];
itaepol1=roots(itaeeq1);
pi1=itaepol1(1);
pi2=itaepol1(2);
pi3=itaepol1(3);
pi4=itaepol1(4);
%Integral ITAE
itaeeq2=[1 2.8*wi 5.0*wi^2 5.5*wi^3 +3.4*wi^4 wi^5];
itaepol2=roots(itaeeq2);
pin1=itaepol2(1);
pin2=itaepol2(2);
pin3=itaepol2(3);
pin4=itaepol2(4);
pin5=itaepol2(5);

%Control implementation
ki=place(A,B,[pi1,pi2,pi3,pi4]);
kpi=1/-((C-D*ki)*inv(A-B*ki)*B-D);
plant_ucitae=ss(A-B*ki, B*kpi, C-D*ki, D*kpi);

%ITAE with static error
%Kpre Gain
% kprei=1/-((C-D*ki)*inv(A-B*ki)*B-D);
% plant_ucitaek=ss(A-B*ki,B*kprei,C-D*ki,D*kprei);
%Integral Part
kiext=place(Aext,Bext,[pin1,pin2,pin3,pin4,pin5]);
ki2=kiext(1:n);
kei=kiext(n+1);
plant_ucitaei=ss([A-B*ki2 -B*kei;-C+D*ki2 D*kei],[zeros(n,1); 1],[C 0],0);

%% State Observer
% Design
PolesObs = p(1:n)*10;
L = place(A',C', PolesObs)'; %ganancia del obs

% Observer Space State
Aobs = A - L*C;
Bobs = [B-L*D L];
Cobs = eye(n);
Dobs = zeros(n,2);

% Observer + Controller Space State
Aoc = A-L*C-B*kn+L*D*kn;
Boc = [B-L*D L];
Coc = -kn;
Doc = [1 0];

% Observer + Controller Space State + Static Control
Aocc = [A-L*C-B*kn+L*D*kn -B*ke+L*D*ke; zeros(1,n) 0];
Bocc = [zeros(n,1) L; 1 -1];
Cocc = [-kn -ke];
Docc = zeros(1,2);

%% table generations 
[y1,t1]=step(plant);plant_r=stepinfo(y1,t1);
[y2,t2]=step(plant_uc);plant_uc_r=stepinfo(y2,t2);
[y4,t4]=step(plant_uci);plant_uci_r=stepinfo(y4,t4);
[y5,t5]=step(plant_ucitae);plant_ucitae_r=stepinfo(y5,t5);
[y7,t7]=step(plant_ucitaei);plant_ucitaei_r=stepinfo(y7,t7);

ControlType={'Plant';'POLE PLACEMENT';'POLE PLACEMENT INTEGRAL';'ITAE';'ITAE INTEGRAL'};
RiseTime={plant_r.RiseTime;plant_uc_r.RiseTime;plant_uci_r.RiseTime;plant_ucitae_r.RiseTime;plant_ucitaei_r.RiseTime};
SettlingTime={plant_r.SettlingTime;plant_uc_r.SettlingTime;plant_uci_r.SettlingTime;plant_ucitae_r.SettlingTime;plant_ucitaei_r.SettlingTime};
Overshoot={plant_r.Overshoot;plant_uc_r.Overshoot;plant_uci_r.Overshoot;plant_ucitae_r.Overshoot;plant_ucitaei_r.Overshoot};

Metrics_from_matlab=table(ControlType,RiseTime,SettlingTime,Overshoot);
disp(Metrics_from_matlab)

%% Simulations
subplot(3,2,[1,2]);
step(plant), grid on, title('PLANT')
subplot(323);
step(plant_uc), grid on, title('PLANT - POLE PLACEMENT')
subplot(324);
step(plant_uci), grid on, title('PLANT - POLE PLACEMENT INTEGRAL')
subplot(325);
step(plant_ucitae), grid on, title('PLANT UNDER CONTROL ITAE')
subplot(326);
step(plant_ucitaei), grid on, title('PLANT - ITAE INTEGRAL')
sgtitle('SIMULATIONS FROM MATLAB')
%% Simulink
%Input and output disturbances
PI=0.000001;
PO=0.000001;
%Parametric uncertainty
INC=25;
A1 = [0 1 0 0;(-c/th1)*INC (-(d1+d)/th1)*INC (-c/th1)*INC (-d/th1)*INC; 0 0 0 1; (-c/th2)*INC (-d/th2)*INC (-c/th2)*INC (-(d2+d)/th2)*INC];

%% run simulink file 
model=sim('Control_Plant_Simulation_2018b')

%% table generations for simulinki 
model_plant_r=stepinfo(model.plant.data,model.plant.time);
model_plant_uc_r=stepinfo(model.plant_uc.data,model.plant_uc.time);
model_plant_r_uc_pu=stepinfo(model.plant_uc_pu.data,model.plant_uc_pu.time);
model_plant_r_uc_iod=stepinfo(model.plant_uc_iod.data,model.plant_uc_iod.time);
model_plant_uci_r=stepinfo(model.plant_uci.data,model.plant_uci.time);
model_plant_ucipu_r=stepinfo(model.plant_uci_pu.data,model.plant_uci_pu.time);
model_plant_uciiod_r=stepinfo(model.plant_uci_iod.data,model.plant_uci_iod.time);
model_plant_ucitae_r=stepinfo(model.plant_ucitae.data,model.plant_ucitae.time);
model_plant_r_ucitae_pu=stepinfo(model.plant_ucitae_pu.data,model.plant_ucitae_pu.time);
model_plant_r_ucitae_iod=stepinfo(model.plant_ucitae_iod.data,model.plant_ucitae_iod.time);
model_plant_ucitaei_r=stepinfo(model.plant_ucitaei.data,model.plant_ucitaei.time);
model_plant_ucitaepu_r=stepinfo(model.plant_ucitae_pu.data,model.plant_ucitae_pu.time);
model_plant_ucitaeiiod_r=stepinfo(model.plant_ucitae_iod.data,model.plant_ucitae_iod.time);

ControlType={'Plant';'POLE PLACEMENT';'POLE PLACEMENT PARAM. UNCERTAINTY';'POLE PLACEMENT I/O DISTURB.';'POLE PLACEMENT INTEGRAL';'POLE PLACEMENT INTEGRAL PARAM. UNCERTAINTY';'POLE PLACEMENT INTEGRAL I/O DISTURB.';'ITAE';'ITAE PARAM. UNCERTAINTY';'ITAE I/O DISTURB.';'ITAE INTEGRAL';'ITAE INTEGRAL PARAM. UNCERTAINTY';'ITAE INTEGRAL I/O DISTURB.'};
RiseTime={model_plant_r.RiseTime;model_plant_uc_r.RiseTime;model_plant_r_uc_pu.RiseTime;model_plant_r_uc_iod.RiseTime;model_plant_uci_r.RiseTime;model_plant_ucipu_r.RiseTime;model_plant_uciiod_r.RiseTime;model_plant_ucitae_r.RiseTime;model_plant_r_ucitae_pu.RiseTime;model_plant_r_ucitae_iod.RiseTime;model_plant_ucitaei_r.RiseTime;model_plant_ucitaepu_r.RiseTime;model_plant_ucitaeiiod_r.RiseTime};
SettlingTime={model_plant_r.SettlingTime;model_plant_uc_r.SettlingTime;model_plant_r_uc_pu.SettlingTime;model_plant_r_uc_iod.SettlingTime;model_plant_uci_r.SettlingTime;model_plant_ucipu_r.SettlingTime;model_plant_uciiod_r.SettlingTime;model_plant_ucitae_r.SettlingTime;model_plant_r_ucitae_pu.SettlingTime;model_plant_r_ucitae_iod.SettlingTime;model_plant_ucitaei_r.SettlingTime;model_plant_ucitaepu_r.SettlingTime;model_plant_ucitaeiiod_r.SettlingTime};
Overshoot={model_plant_r.Overshoot;model_plant_uc_r.Overshoot;model_plant_r_uc_pu.Overshoot;model_plant_r_uc_iod.Overshoot;model_plant_uci_r.Overshoot;model_plant_ucipu_r.Overshoot;model_plant_uciiod_r.Overshoot;model_plant_ucitae_r.Overshoot;model_plant_r_ucitae_pu.Overshoot;model_plant_r_ucitae_iod.Overshoot;model_plant_ucitaei_r.Overshoot;model_plant_ucitaepu_r.Overshoot;model_plant_ucitaeiiod_r.Overshoot};

Metrics_from_simulink=table(ControlType,RiseTime,SettlingTime,Overshoot);
disp(Metrics_from_simulink)

%% Simulations
subplot(4,4,[1,2,3,4]);
plot(model.plant), grid on, title('PLANT')
subplot(445);
plot(model.plant_uc), grid on, title('PLANT - POLE PLACEMENT')
subplot(446);
plot(model.plant_uc_pu), grid on, title('POLE PLACEMENT PARAM. UNCERTAINTY')
subplot(447);
plot(model.plant_uc_iod), grid on, title('POLE PLACEMENT I/O DISTURB.')
subplot(448);
plot(model.plant_uci), grid on, title('PLANT - POLE PLACEMENT INTEGRAL')
subplot(449);
plot(model.plant_uci_pu), grid on, title('PLANT - POLE PLACEMENT INTEGRAL PARAM. UNCERTAINTY')
subplot(4,4,10);
plot(model.plant_uci_iod), grid on, title('PLANT - POLE PLACEMENT INTEGRAL I/O DISTURB.')
subplot(4,4,11);
plot(model.plant_ucitae), grid on, title('PLANT - ITAE')
subplot(4,4,12);
plot(model.plant_ucitae_pu), grid on, title('PLANT - ITAE PARAM. UNCERTAINTY ')
subplot(4,4,13);
plot(model.plant_ucitae_iod), grid on, title('PLANT - ITAE PARAM. I/O DISTURB. ')
subplot(4,4,14);
plot(model.plant_ucitaei), grid on, title('PLANT - ITAE INTEGRAL')
subplot(4,4,15);
plot(model.plant_ucitae_pu), grid on, title('PLANT - ITAE INTEGRAL PARAM. UNCERTAINTY ')
subplot(4,4,16);
plot(model.plant_ucitae_iod), grid on, title('PLANT - ITAE INTEGRAL I/O DISTURB.')

sgtitle('SIMULATIONS FROM SIMULINK')
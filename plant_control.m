%% Definition of the Plant
%Variables and constants:
th1=0.0019;
th2=0.0019;
kt1=3.6242*exp(-005);
kt2=2.05*exp(-005);
c=0.2038;
d=0.0022;
d1=4.7*exp(005);
d2=4.7*exp(005);

%Plant:
A=[0 1 0 0;-c/th1 -(d1+d)/th1 -c/th1 -d/th1;0 0 0 1; -c/th2 -d/th2 -c/th2 -(d2+d)/th2];
B=[0 0;kt1/th1 0; 0 0; 0 kt2/th2];
C=[1 0 0 0; 0 0 1 0];
D=[0 0;0 0];
Plant=ss(A,B,C,D);
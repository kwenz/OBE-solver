clear;

import bloch.*

warning('off','all');


syms wg w0 we D real;
syms W_L w_L W_m w_m real;
syms d_L d_m real;
syms g real;
assume(g,'positive')






H=Hamiltonian(4);
H.addEnergies([w0,wg,wg+D,we]);
H.addCoupling(1,2,W_m,w_m);
H.addCoupling(1,3,W_m,w_m);
H.addCoupling(2,4,W_L,w_L);
H.addCoupling(3,4,W_L,w_L);
H.defineZero(we);
H.unitaryTransformation();
H.defineEnergyDetuning(w0,wg,d_m,w_m);
H.defineEnergyDetuning(wg,we,d_L,w_L);

F=H.transformed;
v1=[1;0;0;0];
v2=[0;1/sqrt(2);1/sqrt(2);0];
v3=[0;1/sqrt(2);-1/sqrt(2);0];
v4=[0;0;0;1];
U=[v1,v2,v3,v4];

F=simplify(U'*F*U,'Steps',100);
% F=subs(F,d_L,D/2);
% F=subs(F,d_m,-D/2);
disp(F)

% G=F(1:2,1:2);
% 
% [A,B]=eig(G);
% 
% % disp(A)
% % disp(B)
% 
syms a1 a2 b1 b2 real;
% 
Ww=sqrt(W_m^2+D^2/2);
% Ww=1;
x1=[W_m/Ww;0;D/(sqrt(2)*Ww);0];
x2=[0;1;0;0];
x3=[D/(sqrt(2)*Ww);0;-W_m/Ww;0];
x4=[0;0;0;1];
X=[x1,x2,x3,x4];

disp(U*X)


F=simplify(X'*F*X,'Steps',100);
% F=subs(F,
disp(F)
return
% H.transformed=F;
L=Dissipator(4);
L.addDecay(4,3,g/2);
L.addDecay(4,2,g/2);


% H.createGraph(L);
% H.plotGraph();
gamma=1;
Omega_m=2;
Omega_L=10;
Delta=10;
del_m=0;
del_L=Delta/2;

eq=BlochEqns(H,L);
IC=zeros(4);
% 
IC(1,1)=0;
IC(2,2)=1/2;
IC(3,3)=1/2;
% 
% eq.necessaryVariables();
% return

eq.evolve(0,100,IC,[Delta,Omega_L,Omega_m,del_L,del_m,gamma]);
eq.plotEvolution();
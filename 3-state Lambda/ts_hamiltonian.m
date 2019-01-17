clear all;
import bloch.*

warning('off','all');


% Hamiltonian

syms wg we g B wx1 wx2 wz Wz Wx1 Wx2 W a D dz dx1 dx2 Bet1 Bet2 ws1 ws2 real;

assume(W,'positive');
H=Hamiltonian(4);
H.addEnergies([wg-g*B,wg,wg+g*B,we]);


% h=H.hamiltonian;
% 
% h(1,3)=h(1,3)+g*B;
% h(3,1)=h(3,1)+g*B;
% 
% H.hamiltonian=h;

H.addCoupling(2,4,Wz,wz);
H.addCoupling(1,4,Wx1,wz);
H.addCoupling(3,4,Wx2,wz);
H.defineZero(we);
H.unitaryTransformation();
H.defineStateDetuning(2,4,dz);
% H.defineStateDetuning(1,4,dx1);
% H.defineStateDetuning(3,4,dx2);
F=H.transformed;




W1=sqrt(Wz^2+Wx1^2);
W2=sqrt(Wx1^2+Wx2^2);
W3=sqrt(Wz^2+Wx2^2);

Ww=sqrt(Wz^2+Wx1^2+Wx2^2);
W0=Ww*W2;


v2=1/W0*[Wx1*Wz;-(Wx1^2+Wx2^2);Wx2*Wz;0];
v1=1/Ww*[Wx1;Wz;Wx2;0];
v3=1/W2*[-Wx2;0;Wx1;0];
v4=[0;0;0;1];

U=[v1,v2,v3,v4];







F=simplify(U'*F*U,'Steps',100);
U=subs(U,[Wx1,Wx2,Wz],[W*sin(a)/sqrt(2),W*sin(a)/sqrt(2),W*cos(a)]);

disp(simplify(U,'Steps',100))

F=subs(F,[Wx1,Wx2,Wz],[W*sin(a)/sqrt(2),W*sin(a)/sqrt(2),W*cos(a)]);
% F=subs(F,[dx1,dx2],[dz+D,dz-D]);

disp(simplify(F,'Steps',200));
return
H.transformed=F;

syms G real;
assume(G,'positive')
d=Dissipator(4);
d.addDecay(4,1,G/3);
d.addDecay(4,2,G/3);
d.addDecay(4,3,G/3);
% 
% 
eq=BlochEqns(H,d);
% 

Om=1;
gfact=1;
B=0;
delx1=2;
delx2=2;
delz=0;
gamma=1;




% eq.necessaryVariables();
% return

IC=zeros(4);
IC(1,1)=1;

eq.evolve(0,200,IC,[gamma,Om,delx1,delx2,delz]);
eq.plotEvolution();

eq.solveSteady();

pee=eq.steadyState(4,4);
disp(simplify(pee,'Steps',500))
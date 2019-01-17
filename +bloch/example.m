clear
import bloch.*

warning('off','all');


%Hamiltonian

syms wa wb we Wae wae Wbe wbe dae dbe real;


h=Hamiltonian(3);
h.addEnergies([wa,wb,we]);
h.addCoupling(1,3,Wae,wae);
h.addCoupling(2,3,Wbe,wbe);
h.defineStateDetuning(1,3,dae);
h.defineStateDetuning(2,3,dbe);
h.defineZero(we);
h.unitaryTransformation();


disp('Hamiltonian before and after unitary transformation')
% disp(h.hamiltonian)
% disp(h.transformed)


%Dissipator

syms Ga Gb real;
assume(Ga,'positive')
assume(Gb,'positive')
d=Dissipator(3);
d.addDecay(3,1,Ga);
d.addDecay(3,2,Gb);

% disp('Dissipator and branching ratios')
% disp(d.dissipator)
% disp(d.decayR)
% disp(d.branching)

%Graph


% h.createGraph(d);
% h.plotGraph();


%Bloch equations

eq=BlochEqns(h,d);


%Steady state solution
% eq.solveSteady();
% pee=eq.steadyState(3,3);
% 
% 
% disp('Steady state solution for the excited state population')
% disp(pee)


% Time evolution

% eq.necessaryVariables(); %Called first, to see which variables need to be substituted.
% return 


IC=zeros(3);  
IC(1,1)=1;

gamma_a=0.6;
gamma_b=0.4;
Rabi_a=1;
Rabi_b=1;
delta_a=0;
delta_b=0;

eq.eqnsRHS=subs(eq.eqnsRHS,[Ga,Gb],[gamma_a,gamma_b]);
%%
Params={Wae,0,1;Wbe,0,1;dae,-5,5;dbe,-5,5};

eq.optimizeParameters(0,50,IC,Params,[3],'Integration','no','Popsize',10,'Iterations',10,'Popnumber',1);
%%
eq.evolve(0,50,IC,eq.optParams);
eq.plotEvolution();




% eq.evolve(0,25,IC,[gamma_a,gamma_b,Rabi_a,Rabi_b,delta_a,delta_b]);
% 
% eq.plotEvolution();




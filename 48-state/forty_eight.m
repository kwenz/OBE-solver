clear;
warning('off','all');
format short;

import bloch.*
import TlF.*

time_of_flight=50; %[us]
switching_freq=1/5;%[MHz]
no_switches=floor(time_of_flight*switching_freq);
t_step=time_of_flight*1000/no_switches; %[ns]

% Physical constants

hbar=1.054*10^(-34); %[Js]
k_b=1.381*10^(-23); % [J/K]
c=299792000; %[m/s]
eps_0=8.854*10^(-12); %[F/m]
Gamma=1/99; %1/Lifetime of the excited state [GHz]
a0=5.29*10^(-11); % [m]
q_e=1.602*10^(-19); % [C]
r_expval=2*pi*0.16*a0; % [m]
B_0=6.686667; % [GHz]
T=10; %[K]
%%
% States
J_g=[0,1,2];
J_e=[1];

St=generateStates(J_g,J_e,0);
n=size(St,1);





% Transition matrix

rabi_matrix=zeros(n,n,3);
disp('Transition matrix')
ind=0;
for p=[-1,0,1]
    ind=ind+1;
    for f=1:n     %f - final states
        for i=1:n %i - initial states
            rabi_matrix(i,f,ind)=dipoleTransitionMatrixElement(St(f,:),St(i,:),p);
        end
    end
end


for i=40:-1:37
    rabi_matrix(i,:,:)=[];
    rabi_matrix(:,i,:)=[];
end

n=n-4;

% Initial distribution
disp('Initial distribution')

Ini=boltzmann_distribution(B_0*10^9,T,0:20);
IC=zeros(n);
ind=0;
for j=J_g    
    for i=1:4*(2*j+1)
        k=ind+i;
        IC(k,k)=Ini(j+1)/(4*(2*j+1));
    end
    ind=ind+4*(2*j+1);
end
IC=IC./sum(Ini(1:max(J_g)+1));



%%

% Variables

syms w_0 w_1 w_2 w_e real; %energies
syms w_m0 w_m1 w_L real; %light frequencies
syms W_m0 W_m1 W_L real; %Rabi frequencies
syms d_m0 d_m1 d_L real; %detunings
syms D_0 D_1 D_2 real; %splittings
syms a_m0 p_m0 a_m1 p_m1 real; % microwave angles: a_m = polarization angle, p_m = beam direction angle with respect to quant axis z 
syms pol_L real; % laser polarization indicator variable: 1 - sigma+, 1 - sigma-

% Their values


E_0=0; % [GHz]
E_1=0;%2*B_0;  
E_2=0;%6*B_0;
E_e=0;%1103400;
Del_e=13.52;
Del_1=0.0004;
Del_2=0.0005;
Del_0=0.0002;
det_L=Del_2/2;
det_m0=1.1*Del_0;
det_m1=0;
polangle_m0=pi/2;
dirangle_m0=pi/4;
polangle_m1=pi/4;
dirangle_m1=-pi/4;

% Light Properties
%%


P_L=0.01; %[W]
diam_L=0.017; %[m]
P_m0=0.00005; %[W]
diam_m0=0.15; %[m]
P_m1=0.005; %[W]
diam_m1=0.1; %[m]

Rabi_m0=10^(-9)*q_e*r_expval*sqrt(8*P_m0/(pi*c*eps_0*diam_m0^2))/hbar; %[GHz]
Rabi_m1=10^(-9)*q_e*r_expval*sqrt(8*P_m1/(pi*c*eps_0*diam_m1^2))/hbar; %[GHz]
Rabi_L=10^(-9)*q_e*r_expval*sqrt(8*P_L/(pi*c*eps_0*diam_L^2))/hbar;

disp(vpa(Rabi_L/Gamma,3))
disp(vpa(Rabi_m0/Gamma,3))
disp(vpa(Rabi_m1/Gamma,3))
%%
% Branching ratios
disp('Branching ratios')


BR=zeros(n);

transition_strengths=zeros(n);
for i=1:n
    for f=1:n
        for p=1:3
            transition_strengths(i,f)=transition_strengths(i,f)+rabi_matrix(f,i,p)^2;
        end
    end
end
for i=1:n
    sums=0;
    for f=1:n
        sums=sums+transition_strengths(i,f);
    end
    for f=1:n
        BR(i,f)=transition_strengths(i,f)/sums;
    end
end

for i=1:n-8
    BR(i,:)=0;
end


% Dissipator

disp('Dissipator')

L=Dissipator(n);

syms G real;
assume(G,'positive') 


DR=zeros(1,n-8);
DR=[DR,G,G,G,G,G,G,G,G];



L.fromBranching(BR,DR);


%%
%Laser
WLp=W_L*(pol_L+1)/2;
WLm=W_L*(1-pol_L)/2;
WLz=0;

WL_pol=[WLm,WLz,WLp];

%Microwaves
Wm0p=W_m0/sqrt(2)*(sin(a_m0)+1i*cos(p_m0)*cos(a_m0));
Wm0m=W_m0/sqrt(2)*(sin(a_m0)-1i*cos(p_m0)*cos(a_m0));
Wm0z=-W_m0*cos(a_m0)*sin(p_m0);

Wm0_pol=[Wm0m,Wm0z,Wm0p];


Wm1p=W_m1/sqrt(2)*(sin(a_m1)+1i*cos(p_m1)*cos(a_m1));
Wm1m=W_m1/sqrt(2)*(sin(a_m1)-1i*cos(p_m1)*cos(a_m1));
Wm1z=-W_m1*cos(a_m1)*sin(p_m1);

Wm1_pol=[Wm1m,Wm1z,Wm1p];

%%
% Hamiltonian

import bloch.*
import TlF.*

disp('Hamiltonian')

H=Hamiltonian(n);

H.addEnergies([w_0,w_0+D_0,w_0+D_0,w_0+D_0,...
    w_1,w_1,w_1,w_1,...
    w_1+D_1,w_1+D_1,w_1+D_1,w_1+D_1,w_1+D_1,w_1+D_1,w_1+D_1,w_1+D_1...
    w_2,w_2,w_2,w_2,w_2,w_2,w_2,w_2,...
    w_2+D_2,w_2+D_2,w_2+D_2,w_2+D_2,w_2+D_2,w_2+D_2,w_2+D_2,w_2+D_2,w_2+D_2,w_2+D_2,w_2+D_2,w_2+D_2,...
    w_e,w_e,w_e,w_e,w_e,w_e,w_e,w_e]);


for i=1:n
    for f=i:n
        for j=1:3
            if rabi_matrix(i,f,j)~=0 && f>n-8 && i<=n-8 && i>4
                Wr=(-1)^(j)*WL_pol(j)*rabi_matrix(i,f,j);
                H.addCoupling(i,f,Wr,w_L);
            elseif rabi_matrix(i,f,j)~=0 && f<=n-8 && i<=n-8 && i>4
                Wr=(-1)^(j)*Wm1_pol(j)*rabi_matrix(i,f,j);
                H.addCoupling(i,f,Wr,w_m1);
            elseif rabi_matrix(i,f,j)~=0 && f<=n-8 && i<=4
                Wr=(-1)^(j)*Wm0_pol(j)*rabi_matrix(i,f,j);
                H.addCoupling(i,f,Wr,w_m0);
            end
        end
    end
end


H.defineEnergyDetuning(w_2,w_e,d_L,w_L);
H.defineEnergyDetuning(w_1,w_2,d_m1,w_m1);
H.defineEnergyDetuning(w_0,w_1,d_m0,w_m0);

H.defineZero(w_2);


H.unitaryTransformation();

disp(H.transformed)


%%
import bloch.*
% Bloch Equations
disp('Optical Bloch Equations')

Eq=BlochEqns(H,L);

Eq.eqnsRHS=subs(Eq.eqnsRHS,[w_0,w_1,w_2,w_e,D_1,D_2,D_0,G,a_m0,p_m0,a_m1,p_m1],[E_0,E_1,E_2,E_e,Del_1,Del_2,Del_0,Gamma,polangle_m0,dirangle_m0,polangle_m1,dirangle_m1]);

disp('Solving')
for sw=1:no_switches
    if mod(sw,2)
        polarization_L=1;
    else
        polarization_L=-1;
    end
    
    

if sw==1
    Eq.evolve(0,t_step,IC,[Rabi_L,Rabi_m0,Rabi_m1,det_L,det_m0,det_m1,polarization_L]);
else
    Eq.intTime=[(sw-2)*t_step,(sw-1)*t_step];
    Eq.lastSol=prevSol;
    Eq.extendEvolution(t_step*sw,[Rabi_L,Rabi_m0,Rabi_m1,det_L,det_m0,det_m1,polarization_L]);
end

prevSol=Eq.lastSol;

fprintf('Done %.2f%% \n',sw/no_switches*100);

% % Eq.plotEvolution();
% figure
% plot(Eq.evTime(1,:),squeeze(Eq.evolution(9,9,:)+Eq.evolution(10,10,:)+Eq.evolution(11,11,:)+Eq.evolution(12,12,:)+Eq.evolution(13,13,:)+Eq.evolution(14,14,:)+Eq.evolution(15,15,:)+Eq.evolution(16,16,:)+Eq.evolution(17,17,:)+Eq.evolution(18,18,:)))
% hold on
% plot(Eq.evTime(1,:),squeeze(Eq.evolution(3,3,:)+Eq.evolution(4,4,:)+Eq.evolution(5,5,:)+Eq.evolution(6,6,:)+Eq.evolution(7,7,:)+Eq.evolution(8,8,:)))
% hold on
% plot(Eq.evTime(1,:),squeeze(Eq.evolution(1,1,:)+Eq.evolution(2,2,:)))
% hold on
% plot(Eq.evTime(1,:),squeeze(Eq.evolution(21,21,:)+Eq.evolution(22,22,:)+Eq.evolution(19,19,:)+Eq.evolution(20,20,:)))
% drawnow

end

% Eq.plotEvolution()

% figure
% plot(Eq.evTime(1,:),squeeze(Eq.evolution(9,9,:)+Eq.evolution(10,10,:)+Eq.evolution(11,11,:)+Eq.evolution(12,12,:)+Eq.evolution(13,13,:)+Eq.evolution(14,14,:)+Eq.evolution(15,15,:)+Eq.evolution(16,16,:)+Eq.evolution(17,17,:)+Eq.evolution(18,18,:)))
% hold on
% plot(Eq.evTime(1,:),squeeze(Eq.evolution(3,3,:)+Eq.evolution(4,4,:)+Eq.evolution(5,5,:)+Eq.evolution(6,6,:)+Eq.evolution(7,7,:)+Eq.evolution(8,8,:)))
% hold on
% plot(Eq.evTime(1,:),squeeze(Eq.evolution(1,1,:)+Eq.evolution(2,2,:)))
% hold on
% plot(Eq.evTime(1,:),squeeze(Eq.evolution(21,21,:)+Eq.evolution(22,22,:)+Eq.evolution(19,19,:)+Eq.evolution(20,20,:)))
% drawnow
%%
figure
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(1,1,:)+Eq.evolution(2,2,:)+Eq.evolution(3,3,:)+Eq.evolution(4,4,:)))
hold on
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(5,5,:)+Eq.evolution(6,6,:)+Eq.evolution(7,7,:)+Eq.evolution(8,8,:)+Eq.evolution(9,9,:)+Eq.evolution(10,10,:)+Eq.evolution(11,11,:)+Eq.evolution(12,12,:)+Eq.evolution(13,13,:)+Eq.evolution(14,14,:)+Eq.evolution(15,15,:)+Eq.evolution(16,16,:)))
hold on
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(17,17,:)+Eq.evolution(18,18,:)+Eq.evolution(19,19,:)+Eq.evolution(20,20,:)+Eq.evolution(21,21,:)+Eq.evolution(22,22,:)+Eq.evolution(23,23,:)+Eq.evolution(24,24,:)+Eq.evolution(25,25,:)+Eq.evolution(26,26,:)+Eq.evolution(27,27,:)+Eq.evolution(28,28,:)+Eq.evolution(29,29,:)+Eq.evolution(30,30,:)+Eq.evolution(31,31,:)+Eq.evolution(32,32,:)+Eq.evolution(33,33,:)+Eq.evolution(34,34,:)+Eq.evolution(35,35,:)+Eq.evolution(36,36,:)))
hold on
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(37,37,:)+Eq.evolution(38,38,:)+Eq.evolution(39,39,:)+Eq.evolution(40,40,:)+Eq.evolution(41,41,:)+Eq.evolution(42,42,:)+Eq.evolution(43,43,:)+Eq.evolution(44,44,:)))
legend('J=0','J=1','J=2','Je=1')
drawnow
%%
figure
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(1,1,:)))
hold on
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(2,2,:)))
hold on
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(3,3,:)))
hold on
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(4,4,:)))
legend('F=0, M_F=0','F=1 M_F=-1','F=1 M_F=0','F=1 M_F=1')
xlabel('Time [us]')
drawnow


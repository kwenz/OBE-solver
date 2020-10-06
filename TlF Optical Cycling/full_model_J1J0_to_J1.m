% clear all;
import bloch.*
import TlF.*

warning('off','all');

time_of_flight=10000; %[ns]
pulse_time=300;%[ns]
switching_freq=1/(2*pulse_time);%[GHz]
no_switches=floor(time_of_flight/pulse_time);
t_step=time_of_flight/no_switches; %[ns]
% t_step=1*1000;

% Physical constants

hbar=1.054*10^(-34); %[Js]
k_b=1.381*10^(-23); % [J/K]
c=299792000; %[m/s]
eps_0=8.854*10^(-12); %[F/m]
Gamma=1/99; %1/Lifetime of the excited state [GHz]
a0=5.29*10^(-11); % [m]
q_e=1.602*10^(-19); % [C]
r_expval=0.32*a0; % [m]
B_0=6.686667; % [GHz]
T=7; %[K]
r_expval_m=5.22*r_expval;

%/////////////////////////////////////////////////////////////////////////%

% States
J_g=[0,1];
J_e=[1];

St=generateStates(J_g,J_e,1);
n=size(St,1);

%/////////////////////////////////////////////////////////////////////////%

% Transition matrix
rabi_matrix=load('TransitionMatrix_EtoF_mixing.mat');
rabi_matrix=rabi_matrix.rabi_matrix;
% 
% rabi_matrix=zeros(n,n,3);
% disp('Transition matrix')
% ind=0;
% for p=[-1,0,1]
%     ind=ind+1;
%     for f=1:n     %f - final states
%         for i=1:n %i - initial states
%             rabi_matrix(i,f,ind)=dipoleTransitionMatrixElement(St(f,:),St(i,:),p,'Mixing','include');
%         end
%     end
% end

%Elimination of states that don't participate in the process
rabi_matrix(end-11,:,:)=[];
rabi_matrix(:,end-11,:)=[];
rabi_matrix(end-7:end,:,:)=[];
rabi_matrix(:,end-7:end,:)=[];
rabi_matrix(17:end-3,:,:)=[];
rabi_matrix(:,17:end-3,:)=[];


n=n-9;

%/////////////////////////////////////////////////////////////////////////%

% Branching ratios
disp('Branching ratios')
BR=load('BranchingRatios_EtoF_mixing.mat');
BR=BR.BR;




%Elimination of excited states that don't participate in the process
BR(end-11,:)=[];
BR(:,end-11)=[];
BR(end-7:end,:)=[];
BR(:,end-7:end)=[];
BR(17:end-3,:)=[];
BR(:,17:end-3)=[];


% BR=zeros(n);
% transition_strengths=zeros(n);
% for i=1:n
%     for f=1:n
%         for p=1:3
%             transition_strengths(i,f)=transition_strengths(i,f)+rabi_matrix(f,i,p)^2;
%         end
%     end
% end
% for i=1:n
%     sums=0;
%     for f=1:n
%         sums=sums+transition_strengths(i,f);
%     end
%     for f=1:n
%         BR(i,f)=transition_strengths(i,f)/sums;
%     end
% end
% 
% for i=1:n-3
%     BR(i,:)=0;
% end

% format rat
% disp(BR)
% format short


%/////////////////////////////////////////////////////////////////////////%

% Dissipator


disp('Dissipator')

L=Dissipator(n);

syms G real;
assume(G,'positive') 


%There are 3 excited states
DR=zeros(1,n-3);
DR=[DR,G,G,G];


L.fromBranching(BR,DR);


%/////////////////////////////////////////////////////////////////////////%

% Variables

import bloch.*

syms w_0 w_1 w_e w_e2 w_3 real; %energies
syms w_m w_L w_L2 real; %light frequencies
syms W_m W_L W_L2 real; %Rabi frequencies
syms d_m d_L d_L2 real; %detunings
syms D_0 D_10 D_1 D_11 real; %splittings
syms a_m a_L real; %polarization angles
syms p_m real; %microwave angle with respect to y (propagation direction of the laser)
syms pol_L pol_L2 pol_m real;

assume(W_L,'positive');
% assume(W_L2,'positive');
assume(W_m,'positive');

%%
%/////////////////////////////////////////////////////////////////////////%

%Hamiltonian
import bloch.*
disp('Hamiltonian')

syms t real;

H=Hamiltonian(n);

H.addEnergies([w_0,w_0-D_0,w_0-D_0,w_0-D_0,...
    w_1,w_1+D_10,w_1+D_10,w_1+D_10,...
    w_1+D_1+D_10,w_1+D_1+D_10,w_1+D_1+D_10,...
    w_1+D_1+D_10+D_11,w_1+D_1+D_10+D_11,w_1+D_1+D_10+D_11,w_1+D_1+D_10+D_11,w_1+D_1+D_10+D_11...
    w_e,w_e,w_e]);

%/////////////////////////////////////////////////////////////////////////%

%Laser
WLp=W_L*(1-sign(sin(2*pi*switching_freq*t)))/2/sqrt(2);
WLm=W_L*(1-sign(sin(2*pi*switching_freq*t)))/2/sqrt(2);
WLz=W_L*(1+sign(sin(2*pi*switching_freq*t)))/2;


WL_pol=[WLm,WLz,-WLp];



Wmp=W_m*exp(1i*pi/4)*(1-sign(sin(2*pi*switching_freq*t+pi/2)))/2/sqrt(2);
Wmm=W_m*exp(-1i*pi/4)*(1-sign(sin(2*pi*switching_freq*t+pi/2)))/2/sqrt(2);
Wmz=W_m*(1+sign(sin(2*pi*switching_freq*t+pi/2)))/2;

Wm0_pol=[Wmm,Wmz,-Wmp];


%/////////////////////////////////////////////////////////////////////////%

%Couplings
for i=1:n
    for f=i:n
        for j=1:3
            if rabi_matrix(i,f,j)~=0 && f>n-3 && i<=n-3
                Wr=(-1)^(j)*WL_pol(4-j)*rabi_matrix(i,f,j);
                H.addCoupling(i,f,Wr,w_L);
            elseif rabi_matrix(i,f,j)~=0 && f<=n-3 && i<=4
                Wr=(-1)^(j)*Wm0_pol(4-j)*rabi_matrix(i,f,j);
                H.addCoupling(i,f,Wr,w_m);
            end
        end
    end
end


H.defineEnergyDetuning(w_1,w_e,d_L,w_L);
H.defineEnergyDetuning(w_0,w_1,d_m,w_m);


H.defineZero(w_0);

H.unitaryTransformation();

disp(H.transformed)

%/////////////////////////////////////////////////////////////////////////%
%%
%Splittings
Del_0=2*pi*13.3*10^(-6);
Del_1=2*pi*176*10^(-6);
Del_11=2*pi*22.24*10^(-6);
Del_12=2*pi*14.54*10^(-6);
%Detunings
det_L=0;
det_m=0;

% Light Properties

P_L=0.01; %[W]
diam_L=0.0012; %[m]
P_m0=1; %[W]
diam_m0=0.03; %[m]

Rabi_m0=10^(-9)*q_e*r_expval_m*sqrt(8*P_m0/(pi*c*eps_0*diam_m0^2))/hbar; %[GHz]
Rabi_L=10^(-9)*q_e*r_expval*sqrt(8*P_L/(pi*c*eps_0*diam_L^2))/hbar;

disp('Rabi rate')
disp('R_L [\Gamma]')
disp(vpa(Rabi_L/Gamma,3))
disp('R_m0 [\Gamma]')
disp(vpa(Rabi_m0/Gamma,3))

%/////////////////////////////////////////////////////////////////////////%

%Initial conditions
IC=zeros(n);

for i=1:4
IC(i,i)=0.0668742;
end
for i=5:16
IC(i,i)=(1-4*0.0668742)/12;
end

% disp(trace(IC))

%/////////////////////////////////////////////////////////////////////////%
%%
%Equations

import bloch.*
Eq=BlochEqns(H,L);

Eq.eqnsRHS=subs(Eq.eqnsRHS,[G,D_0,D_10,D_1,D_11],[Gamma,Del_0,Del_11,Del_1,Del_12]);

Eq.necessaryVariables();

return
%/////////////////////////////////////////////////////////////////////////%
%%
%Solutions

fprintf('R_L = %.2f \n',R_L);
fprintf('R_m = %.2f \n',R_m);
Rabi_L=R_L*Gamma;
Rabi_m=R_m*Gamma;

disp('Solving')
Eq.evolve(0,time_of_flight,IC,[Rabi_L,Rabi_m,0,0]);
disp('Solved')
        
%Excited state evolution
pee=real(squeeze(Eq.evolution(17,17,:))+squeeze(Eq.evolution(18,18,:))+squeeze(Eq.evolution(19,19,:)));

%Total signal, i.e. integral of excited state population over interaction
%time
total_int=trapz(Eq.evTime(:),pee(:));% [ns]

%Average scattering rate
rate=total_int*Gamma/(Eq.evTime(end)); 

%Average scattering rate in units of Gamma
rate_gamma=total_int/(Eq.evTime(end)); 

%Denominator of scattering rate in units of Gamma, e.g.
%rate_gamma=Gamma/10, invrate_gamma=10
invrate_gamma=Gamma/rate;

%Total photon number
ph_num=total_int*Gamma;

fprintf('Int_pee = %.2f \n',total_int);
fprintf('R_avg = %.2f GHz\n',rate);
fprintf('R_avg = %.2f \Gamma\n',rate_gamma);
fprintf('R_avg denominator = %.2f \n',invrate_gamma);
fprintf('N_ph = %.2f \n',ph_num);


%All excited states only plot
figure
plot(Eq.evTime(1,:)./1000,pee)
xlabel('Time [\mu s]')
ylabel('Excited State Population')
drawnow
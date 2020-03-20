% clear all;
import bloch.*
import BaH.* %Small module that includes state generation, dipole transition matrix elements calculation and Zeeman effect calculation in Barium Hydride

warning('off','all');


%***********************PHYSICAL CONSTANTS********************************%

%Frequency units used in this calculation are GHz, and time units of ns.

time_of_flight=20; %[us]

hbar=1.054*10^(-34); %[Js]
k_b=1.381*10^(-23); % [J/K]
c=299792000; %[m/s]
eps_0=8.854*10^(-12); %[F/m]
Gamma=0.964/136.5; %1/Lifetime of the excited state [GHz]
a0=5.29*10^(-11); % [m]
q_e=1.602*10^(-19); % [C]
r_expval=1.632*a0; % [m]
B_field=2*pi*9; %[G]
g_factor=-0.5/1000; %[GHz/G]


%*************************STATE GENERATION********************************%

[StG,StE]=generateStates(); %For reference, check +BaH/generateStates.m

n=12+4; %Total number of states


%*************************TRANSITION DIPOLES******************************%

disp('Transition matrix')

rabi_matrix=zeros(n,n,3);

ind=0;
for p=[-1,0,1] %Three basic polarizations: -1 - right circular, +1 - left circular, 0 - pi
    ind=ind+1;
    for f=1:4     %f - excited electronic states
        for i=1:12 %i - ground electronic states
            rabi_matrix(i,12+f,ind)=dipoleTransitionMatrixElement(StG(i,:),StE(f,:),p); %For reference, +BaH/dipoleTransitionMatrixElement.m
            rabi_matrix(12+f,i,ind)=dipoleTransitionMatrixElement(StE(f,:),StG(i,:),p);            
        end
    end
end


%*************************ZEEMAN EFFECT***********************************%

disp('Zeeman effect matrix')

zeeman_matrix=zeros(n,n);

syms g_f B real; %Symbolic variables for g-factor and magnetic field


%We calculate the Zeeman effect directly only for the ground state, which
%is in Hund's case (b) and for which we found appropriate formulas
for i=1:12     %f - final states
   for f=1:12 %i - initial states
       zeeman_matrix(i,f)=zeemanElement(StG(i,:),StG(f,:)); %For reference, +BaH/zeemanElement.m
    end
end


%*************************BRANCHING RATIOS********************************%


disp('Branching ratios')

BR=zeros(n);

%To calculate branching ratios, we first add squares of all transition dipoles
%(so transition strengths), and then divide appropriate transition
%strengths by that sum.

%Transition strengths
transition_strengths=zeros(n);
for i=1:n
    for f=1:n
        for p=1:3
            transition_strengths(i,f)=transition_strengths(i,f)+rabi_matrix(i,f,p)^2;
        end
    end
end

%Sums of transition strengths for a given initial state 'i'
for i=1:n
    sums=0;
    for f=1:n
        sums=sums+transition_strengths(i,f);
    end
    for f=1:n
        BR(i,f)=transition_strengths(i,f)/sums; %(rotational) branching ratio
    end
end

%Initial states don't decay, so we remove those terms. Otherwise, they
%would simply indicate fractional transition strengths.
for i=1:n-4
    BR(i,:)=0;
end

disp(BR(n-3:end,:)) %Rows indicate branching ratios of decays from a single excited state (columns are ground states)

%****************************DISSIPATOR***********************************%

disp('Dissipator')

L=Dissipator(n);

syms G real;
assume(G,'positive') 

%All excited states are assumed to have the same total decay rate G.
DR=zeros(1,n-4);
DR=[DR,G,G,G,G];

%We use branching ratios table to generate the whole dissipator, instead of
%adding decays one-by-one. That's also why we needed DR vector.
L.fromBranching(BR,DR);


%****************************HAMILTONIAN**********************************%

%Symbolic variables
syms w_0 w_1 w_2 w_e real; %energies
syms w_L01 w_L1 w_L2 real; %light frequencies
syms W_L01 W_L1 W_L2 real; %Rabi frequencies.
syms d_L01 d_L1 d_L2 real; %detunings
syms D_01 D_e real; %splittings
syms a_L01 a_L1 a_L2 real; %polarization angles


%We assume Rabi rates to be real, because the only are related to the 
%electric field magnitude. The complex part originates in polarization 
%vector and is taken into account when adding appropriate couplings to the
%hamiltonian.
assume(W_L01,'positive');
assume(W_L2,'positive');
assume(W_L1,'positive');


%Hamiltonian

disp('Hamiltonian')

H=Hamiltonian(n);


%States are added in following order:
%XSigma
%1 - J=1/2, F=0, m=0
%2-4 - J=1/2, F=1, m=-1,0,1
%5-7 - J=3/2, F=1, m=-1,0,1
%8-12 - J=3/2, F=2, m=-2,-1,0,1,2
%APi
%13 - J=0, F=0, m=0
%14-16 - J=1/2, F=1, m=-1,0,1
H.addEnergies([w_0,w_0-D_01,w_0-D_01,w_0-D_01,...        
    w_1,w_1,w_1,...
    w_2,w_2,w_2,w_2,w_2...
    w_e,w_e+D_e,w_e+D_e,w_e+D_e]);


%Laser 01 - it couples J=1/2 to the excited state. Electric field is
%assumed to have linear polarization at an angle 'a_L01' with respect to
%the quantization axis defined here by the magnetic field. The 'x'
%component decomposes into left and right-circular polarization equally. 
WLp01=W_L01*sin(a_L01)/sqrt(2);
WLm01=W_L01*sin(a_L01)/sqrt(2);
WLz01=W_L01*cos(a_L01);

WL_pol_01=[WLm01,WLz01,-WLp01];


%Laser 1 - it couples J=3/2 F=1 to the excited state. Electric field is
%assumed to have linear polarization at an angle 'a_L1' with respect to
%the quantization axis defined here by the magnetic field.
WLp1=W_L1*sin(a_L1)/sqrt(2);
WLm1=W_L1*sin(a_L1)/sqrt(2);
WLz1=W_L1*cos(a_L1);

WL_pol_1=[WLm1,WLz1,-WLp1];


%Laser 2 - it couples J=3/2 F=2 to the excited state. Electric field is
%assumed to have linear polarization at an angle 'a_L2' with respect to
%the quantization axis defined here by the magnetic field.
WLp2=W_L2*sin(a_L2)/sqrt(2);
WLm2=W_L2*sin(a_L2)/sqrt(2);
WLz2=W_L2*cos(a_L2);

WL_pol_2=[WLm2,WLz2,-WLp2];



%Couplings. 'The W*_pol_*' vectors represent light polarization components.
%The true Rabi rate that appears in the off-diagonal cells that represent
%light coupling, are calculated using a scalar product of rank-1 tensors -
%-T(d)*T(E). This, can be written as: \Sum_{p=[-1,0,1]} (-1)^p * E_{-p} * <f|T(d)_p|i>
for i=1:n
    for f=i:n
        for j=1:3
            %XSigma J=1/2 to APi J=1/2
            if rabi_matrix(i,f,j)~=0 && f>n-4 && i<=4 
                Wr=(-1)^(j)*WL_pol_01(4-j)*rabi_matrix(i,f,j);
                H.addCoupling(i,f,Wr,w_L01);
             
            %XSigma J=3/2 F=1 to APi J=1/2
            elseif rabi_matrix(i,f,j)~=0 && f>n-4 && i<=7 && i>=5
                Wr=(-1)^(j)*WL_pol_1(4-j)*rabi_matrix(i,f,j);
                H.addCoupling(i,f,Wr,w_L1);
                
            %XSigma J=3/2 F=2 to APi J=1/2
            elseif rabi_matrix(i,f,j)~=0 && f>n-4 && i<=12 && i>=8
                Wr=(-1)^(j)*WL_pol_2(4-j)*rabi_matrix(i,f,j);
                H.addCoupling(i,f,Wr,w_L2);
            end
        end
    end
end


%Detunings for all transitions
H.defineEnergyDetuning(w_0,w_e,d_L01,w_L01);
H.defineEnergyDetuning(w_1,w_e,d_L1,w_L1);
H.defineEnergyDetuning(w_2,w_e,d_L2,w_L2);

%Zero is defined as the energy of the excited state
H.defineZero(w_e);


%SIMPLIFACTION USED FOR OPTIMAZTION

% %Unitary Transformation eliminating time dependence
% H.unitaryTransformation();
% % 
% %Next, we add the effects of the magnetic field. Normally, we should add
% %them to the hamiltonian before any transformation is made. In this system,
% %that would result in quickly oscillating (at GHz frequencies) terms
% %proportional to differneces between laser frequencies, e.g. proportional 
% %to g_f*B*exp[i*t*(w_L01-w_L1)], which for w_L01-w_L1 and w_L01-w_L2 is 
% %equal to ~1200 Gamma. The w_L1-w_L2 term is around 40MHz, which is equal
% %to 5.6 Gamma. 
% %
% %The very-quickly oscillating term does not influence our results -
% %oscillations are on top of the much slower behaviour. The ~40MHz
% %oscillating part does affect the transient behaviour in the first few or 
% %tens of microseconds, but at the steady state it only adds oscillations of
% %the same order on top of the behaviour found without these terms. Overall, the
% %average scattering rate over considered interaction times is mostly
% %unaffacted. Lack of those oscillatory terms allows for quicker numeical
% %solutions, and therefore for quicker optimization.
% H.transformed=H.transformed+B*g_f*zeeman_matrix; %Ground states



%PROPER FULL SOLUTION
%In order to obtain a proper hamiltonian, comment out the previous two
%code and uncomment the lines below.

%Zeeman part is added to the hamiltonian in the original basis
H.hamiltonian=H.hamiltonian+B*g_f*zeeman_matrix;

%Unitary transformation is performed
H.unitaryTransformation();

%Next, we change variables to include the real spin-rotation splitting
%between J=1/2 and J=3/2 (D1 = 8.6 GHz) and real splitting between F=0 and
%F=1 states in J=3/2 manifold (D12 = 40 MHz). We do this, because initially
%all three state manfiolds were assumed to have their own energies w_0, w_1
%and w_2. After the transformation terms w_2-w_0, w_1-w_0 and w_2-w_1
%appear in the magnetic couplings (the mentioned quickly oscillating
%terms).
syms D1 D12 real;

H.transformed=simplify(subs(H.transformed,[w_1,w_2],[w_0+D1,w_0+D1+D12]));

%Next, we simplify hamiltonian again. MATLAB sometimes struggles with it.
for i=1:n
    for j=1:n
        if i~=j
%             H.transformed(i,j)=simplify(expand(H.transformed(i,j))); %Use this line if you want to include the fastest oscillating terms (~8.6 GHz)
            H.transformed(i,j)=subs(simplify(expand(H.transformed(i,j))),D1,0); %Use this line to eliminate the fastest oscillating terms.
        end        
    end
end




%Excited states are assumed to be simply shifted in energy. In other words,
%while in the ground states we keep the coupled basis that normally appears 
%without the magnetic field (to see the effects of remixing), for the 
%excited states we assume these are the eigenstates, and so their energy is
%just shifted by B*g_factor. We also assume all g_factors are the same,
%which is true to a good approximation.
H.transformed(14,14)=H.transformed(14,14)+B*g_f; 
H.transformed(16,16)=H.transformed(16,16)-B*g_f;
% 
disp(H.transformed)


%**************************NUMERICAL VALUES*******************************%

%Splittings. Everything is in units of 2\pi*frequency
Del_01=2*pi*2*10^(-3); 
Del_e=2*pi*2*10^(-3);
Del_12=2*pi*40*10^(-3);
Del_1=2*pi*8.6;

%Detunings
det_L01=0.49*Gamma;
det_L1=-5*Gamma;
det_L2=-5*Gamma;

%Polarization angles
angle_L01=1.021;
angle_L1=angle_L01;
angle_L2=angle_L1+pi/2;

%Laser light 
P_L01=0.0511; %[W]
diam_L01=0.01; %[m]
P_L1=0.109; %[W]
diam_L1=0.01; %[m]
P_L2=0.109; %[W]
diam_L2=0.01; %[m]

%Average Rabi rates [GHz]
Rabi_L01=10^(-9)*q_e*r_expval*sqrt(8*P_L01/(pi*c*eps_0*diam_L01^2))/hbar;
Rabi_L1=10^(-9)*q_e*r_expval*sqrt(8*P_L1/(pi*c*eps_0*diam_L1^2))/hbar;
Rabi_L2=10^(-9)*q_e*r_expval*sqrt(8*P_L2/(pi*c*eps_0*diam_L2^2))/hbar;

disp('Rabi rates')
disp('R_L0 [\Gamma]')
disp(vpa(Rabi_L01/Gamma,5))
disp('R_L1 [\Gamma]')
disp(vpa(Rabi_L1/Gamma,5))
disp('R_L2 [\Gamma]')
disp(vpa(Rabi_L2/Gamma,5))


%*************************INITIAL CONDITIONS******************************%

IC=zeros(n);

%We simply assume that 1/4 of population is in J=1/2 and 3/4 in J=3/2. Both
%equally distributed among the Zeeman sublevels
for i=1:4
    IC(i,i)=1/16;
end
for i=5:12
    IC(i,i)=3/32;
end


%*************************MASTER EQUATION*********************************%

disp('Optical Bloch Equations')
Eq=BlochEqns(H,L);

%We substitute values for contant parameters that won't be optimized
Eq.eqnsRHS=subs(Eq.eqnsRHS,[G,D_01,D_e,g_f,D1,D12],[Gamma,Del_01,Del_e,g_factor,Del_1,Del_12]);
Eq.necessaryVariables();


%Example solution
disp('Solving')
Eq.evolve(0,time_of_flight*500,IC,[B_field,Rabi_L01,Rabi_L1,Rabi_L2,angle_L01,angle_L1,angle_L2,det_L01,det_L1,det_L2]);
disp('Solved')
Eq.plotEvolution()

%Excited state evolution
pee=real(squeeze(Eq.evolution(13,13,:))+squeeze(Eq.evolution(14,14,:))+squeeze(Eq.evolution(15,15,:))+squeeze(Eq.evolution(16,16,:)));

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
fprintf('R_avg = %.5f GHz\n',rate);
fprintf('R_avg = %.2f Gamma\n',rate_gamma);
fprintf('R_avg denominator = %.2f \n',invrate_gamma);
fprintf('N_ph = %.2f \n',ph_num);

%All excited states only plot
figure
plot(Eq.evTime(1,:)./1000,pee)
xlabel('Time [\mu s]')
ylabel('Excited State Population')
drawnow
return

%%
%*************************OPTIMIZATION************************************%
import bloch.*
%This section was used for optimization of the scattering rate given the
%experimental constraints.

%We constrain the Rabi rates (W_L2=W_L1), detunings (d_L1=d_L2) and 
%polarization angles (a_L01=a_L1, a_L2=a_L1+\pi/2)
Eq.eqnsRHS=subs(Eq.eqnsRHS,[W_L2,a_L01,a_L2,d_L2],[W_L1,a_L1,a_L1+pi/2,d_L1]);
Eq.necessaryVariables();


%Optimization with respect to: magnetic field (multiplied by 2\pi to obtain
%correct 2\pi*GHz units, Rabi rates, detunings and one free polarization 
%angle. Done over 20us. Optimization performed tomaximize the integral of 
%excited state populations 13-16 over the interaction time.
Params={B,0,2*pi*10;W_L01,0,20*Gamma;W_L1,0,20*Gamma;a_L1,-pi/2,pi/2;d_L01,-10*Gamma,10*Gamma;d_L1,-10*Gamma,10*Gamma};
Eq.optimizeParameters(0,time_of_flight*500,IC,Params,[13,14,15,16],'Iterations',12,'Integration','yes');


%Optimal parameters
fprintf('B = %.2f G\n',Eq.optParams(1)/(2*pi));
fprintf('W_L01 = %.2f Gamma\n',Eq.optParams(2)/Gamma);
fprintf('W_L1 = %.2f Gamma\n',Eq.optParams(3)/Gamma);
fprintf('d_L01 = %.2f Gamma\n',Eq.optParams(5)/Gamma);
fprintf('d_L1 = %.2f Gamma\n',Eq.optParams(6)/Gamma);
fprintf('a_L1 = %.2f pi rad\n',Eq.optParams(4)/pi);


%Time evolution using optimal parameters
Eq.evolve(0,time_of_flight*1000,IC,Eq.optParams);
Eq.plotEvolution()


%Excited state evolution
pee=real(squeeze(Eq.evolution(13,13,:))+squeeze(Eq.evolution(14,14,:))+squeeze(Eq.evolution(15,15,:))+squeeze(Eq.evolution(16,16,:)));

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



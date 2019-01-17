clear;


% Symbolic variables
syms w_a w_b w_e phi real; %Level energies
syms w_ae w_be real; %Light frequencies
syms W_ae W_be real; %Rabi frequencies
syms D_ae D_be d real; %Detunings
t=sym('t','real');
% w_a=-d_ae+w_e;

w_ae=-(D_ae-w_e+w_a); %D_ae=w_ae -(w_a-w_e)  %D_ae=w_e-w_a-w_ae
w_be=-(D_be-w_e+w_b);
% w_b=w_a+d;
% w_be=w_ae ;%Hyperfine splitting, one laser!


%Numerical values of the variables
d_en=2;
E_a=-300; %Energies
E_b=-300+d_en;
E_e=0; %Set to 0

Gamma_a=0.5; %Decay rates
Gamma_b=0.5;

Delta_ae=0; 
Delta_be=0;

Rabi_ae=5; %Rabi frequencies
Rabi_be=5;



%Hamiltonian. We create it first as matrix of symbolic variables, and then
%populate appropriate elements.
H=sym('H',3);
H(:,:)=sym(0);
H(1,1)=w_a;
H(1,3)=W_ae/2*exp(+i*w_ae*t);
H(3,1)=W_ae/2*exp(-i*w_ae*t);
H(2,2)=w_b;
H(2,3)=W_be/2*exp(+i*w_be*t);
H(3,2)=W_be/2*exp(-i*w_be*t);
H(3,3)=w_e;

%Unitary transformation matrix. Diagonal terms can determined by looking 
%at the transformed hamiltonian and demanding simple form for off-diagonal
%elements. 
T=sym('T',3);
T(:,:)=sym(0);
T(1,1)=exp(i*w_ae*t);
T(2,2)=exp(i*w_be*t);
T(3,3)=1;


%Unitary transformation of the hamiltonian.
H_f=T'*H*T-1i*T'*diff(T,t);
H_f=simplify(H_f);

disp(H_f)




% % ------------------------------------------------------------------------------
% % ODE - full time dependent solution


%First, we create symbolic functions representing the density matrix
syms paa(t) pbb(t) pee(t)
syms pab(t) pae(t) 
syms pbe(t) pba(t)
syms peb(t) pea(t)
syms Ga Gb real %Detunings
assume(Ga,'positive')
assume(Gb,'positive')

%Density matrix
P=[paa(t), pab(t), pae(t); pba(t), pbb(t), pbe(t); pea(t), peb(t), pee(t)];

%Density matrix is hermitian 
for i=1:3
    for j=i+1:3
        P(i,j)=conj(P(j,i));
    end
end

%The dissipation matrix
L=sym('L',3);
L(:,:)=sym(0);
L(3,3)=-(Ga+Gb)*pee(t); %Spontaneous decay
L(2,2)=Gb*pee(t);
L(1,1)=Ga*pee(t);

%Off-diagonal elements
GGa(:,:)=sym(0);
GGb(:,:)=sym(0);
GGa(3,3)=sqrt(Ga);
GGb(3,3)=sqrt(Gb);

L=L-1/2*anticommute(GGa*GGa,P)+GGa*P*GGa-1/2*anticommute(GGb*GGb,P)+GGb*P*GGb;


%Now, we create a final matrix, being the right side of Bloch equations...
% Final_matrix=simplify(-1i*commute(H_f,P)+L,'Steps',50);
% 
% %...and it into a vector, where first 3 elements belong to the first row
% %etc.
% Equations=to_vector(Final_matrix);
% n=size(Equations);
% n=n(1);
% 
% 
% 
% %We do the same with the density matrix.
% Variables=to_vector(P);

V1=[w_a,w_b,w_e,D_be,D_ae,W_ae,W_be,Ga,Gb];
V2=[E_a,E_b,E_e,Delta_be,Delta_ae,Rabi_ae,Rabi_be,Gamma_a,Gamma_b];

%Initial conditions
InitialConditions=[1;0;0;0;0;0;0;0;0];

%Time of integration
t_start=0;
t_end=50;

%In case we want to find how the solution depends on the detuning, we loop
%over its values. 
% Doublet=[]; %Results
% tic

ns=brownian_motion(1001,t_end,0,1);
ns=ns(:,1);
% hcn = dsp.ColoredNoise('InverseFrequencyPower',1,'SamplesPerFrame',1001);
% ns = 1/sqrt(1001/t_end)*step(hcn);
figure
plot(0:t_end/1000:t_end,ns);
% return


%     Y=sym('Y', [n 1]); %Symbolic variables used by the ode solver
%     
%     %We substitute symbolic variables with known values and symbolic
%     %variables used by the solver
%     F=subs(Equations,[Variables,w_a,w_b,w_e,D_be,D_ae,W_ae,W_be,Ga,Gb],[Y,E_a,E_b,E_e,Delta_be,Delta_ae,Rabi_ae,Rabi_be,Gamma_a,Gamma_b]);
%     
    %Matlab requires to change the equations in vector form into a matlab
    %function handle
%     FF=matlabFunction(F,'vars',{t,Y});
    
    
    %Options for the solver:
    option=odeset('AbsTol',1e-6,'RelTol',1e-3,'MaxStep',(t_end-t_start)/1000);

    %Solution
    solution=ode45(@(t,Y) defODE(ns,t_end,t,H_f,L,P,V1,V2,Y),[t_start t_end],InitialConditions,option);
    figure
    plot(solution.x(1,:),solution.y(9,:))
   
%     drawnow
%     hold on
%     Doublet=[Doublet,[Delta_ae;solution.y(3,end)]];


% toc

%Matrix to vector
function[vec]=to_vector(A)
    k=size(A,1);
    vec=sym('v',[k^2 1]);
    for i=1:k
        vec(k*(i-1)+1:k*i)=A(:,i);
    end
end

function[X]=commute(A,B)
    X=A*B-B*A;
end

function[X]=anticommute(A,B)
    X=A*B+B*A;
end


function[F]=defODE(noise,tf,t,H_f,L,P,V1,V2,Y)
    T=0:tf/1000:tf;
    disp(t/tf*100)
    ran_phi=interp1(T,noise,t);
    
    H_f(1,3)=H_f(1,3).*exp(1i.*ran_phi);
    H_f(3,1)=conj(H_f(1,3));
    
    Final_matrix=simplify(-1i*commute(H_f,P)+L,'Steps',50);


    Equations=to_vector(Final_matrix);
    n=size(Equations);
    n=n(1);

    Variables=to_vector(P);
        F=subs(Equations,V1,V2);
    
    for i=1:9
        F=subs(F,Variables(i),Y(i));
    end
    
    
    F=double(F);

end


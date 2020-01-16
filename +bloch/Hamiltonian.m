%
% This class defines system's hamiltonian using symbolic variables. If numerical values are necessary, one 
% can always use function subs(obj.***,symbolic variables,numerical values) and store the result 
% in a new user-defined variable. To find symbolic variables that are being used, use symvar(obj.***).
%   
%
%
%
% ///////////////////////////////////////////////////////////////////////////////////////
%
% The stored public properties are:
%
% -> hamiltonian - containes the system's basic hamiltonian. If the hamiltonian is transformed
% using a unitary transformation, it will be stored in a different variable. This one always
% contains the original hamiltonian
%
% -> transformed - containes the hamiltonian after a unitary transformation is applied. For
% basic couplings between states (*exp(iwt)), it doesn't contain their time dependence.
%
% -> transMatrix - unitary transformation matrix used
%
% -> energies - diagonal elements of hamiltonian (energies of all the states) stored as a vector
%
% -> couplings - a NxNx2 matrix (where N is number of states), where (i,j,2) element is Rabi rate associated
% with the i->j transition, and (i,j,1) is frequency of the light field
%
% -> detunings - Mx3 matrix, where M is number of defined detunings. (m,1) is the initial state index, (m,2) the final state index 
% and (m,3) is the detuning d
% 
% -> stateGraph - a directed graph of the system, where edges go along defined couplings (if coupling is defined as i->f,
% the edge will follow that direction) and along decay paths
%
% -> zeroEnergy - energy level (symbolic) defined as having 0 energy. Has to be a part of one of the energy levels specified in 'energies'.
%
%
% The stored private properties are:
%
% -> freqs - stores a vector containing all light frequencies defined
%
% -> cpl - N-element vector, where cpl(i) stores number of couplings from state i
%
% -> cplGr - undirected graph, where edges are couplings between states, which play role of nodes
%
%
% /////////////////////////////////////////////////////////////////////////////////////// 
%
% The class contains following public methods: 
%
% -> contructor Hamiltonian(n) - to create the object in this class, you have to call this function first. 
% The input is number of states (numerical value). The function creates empty matrices and stores them in appropriate variables.
%
% -> function addEnergies(v) - this function takes vector of energies (they should be symbolic) as its input and places them in the 
% hamiltonian. It also stores the vector in the property 'energies'.
%
% -> function addCoupling(S_i,S_f,W_r,w) - couples initial state S_i with final states S_f with a light field
% that has Rabi rate W_r and frequency w. It places appropriate term in the off-diagonal cells of the hamiltonian assuming
% rotating wave approximation (the term is -W_r/2*exp(iwt), without its c.c.; to obtain appropriate Rabi rates with appropriate signs and coefficients
% i.e. coupling strengths, substitute W_r with desired correct value). Variables S_i and S_f are numeric, while W_r and w have to be symbolic!
%
% -> function defineStateDetuning(S_i,S_f,d) - often it is easier to work with detunings than with light frequencies. 
% This function changes variables, by replacing light's frequency w coupling states S_i and S_f with detuning d.
% Variables S_i and S_f are numeric, while d has to be symbolic!
%
% -> function defineEnergyDetuning(E_i,E_f,d,w) - This function changes variables in a more straight-forward way.
% It replaces frequency w by energy difference minus detuning (w=E_f-E-i-d). All variables have to be symbolic.
%
% -> function unitaryTransformation('Method',{'general','equations'}) - The function finds the appropriate frame, where there would be no time dependence in the hamiltonian,
% by using one of two methods. The default one (called without any arguments) looks at coupling through light fields between a state and state
% that is defined as having zero energy. The sum of frequencies in these couplings is put into the unitary matrix and used to transform the Hamiltonian.
% Time dependence is usually eliminated (it is not always possible; it depends on the system) and substitutions from defined detunings are made.
% If the function is called with option ('Method','equation'), the unitary transformation is found by solving a system of linear equation. This method is usable
% if there are at most as many couplings as number states. It does not require a defined zero energy.
% The resulting hamiltonian is stored in property 'transformed' and transformation matrix in 'transMatrix'.
%
% -> function createGraph(L) - uses the dissipator object to create a directed graph. Vertices are states labeled using their respective energies,
% while edges are either couplings between states or decay paths labeld by either Rabi rates or decay rates.
%
% -> function plotGraph() - plots obj.stateGraph of the system
%
% -> function addSidebands(center,order,depth,spacing) - adds 'order' number of sidebands around 'center' frequency. They are added using
% expansion to Bessel functions (so it is crucial that the 'order' is chosen correctly depending on 'depth'; for exact solution, use addPhaseModulation()).
% The modulation depth and spacing between sidebands are define in 'depth' and 'spacing' respectively. This function should be used only after the
% unitary transformation of the hamiltonian (so after unitaryTransformation() function is are called).
%
% -> function addPhaseModulation(center,depth,mod_freq) - adds phase modulation to the laser of 'center' frequency of depth 'depth' and
% modulation frequency 'mod_freq'. The modulating part goes from exp(iwt) to exp(iwt + iBsin(w_m*t)), where 'B' is modulation depth, and 'w_m' is
% modulation frequency.
% 
% **************************************************************************************************************************************************************
% 
% EXAMPLE - 3-level Lambda system
% 
%                                          |e>
%                                     ------------    State 3
%                                  dae ----  ---- dbe
%                                        ^    ^ 
%                                       /      \
%                             Wae, wae /        \ Wbe, wbe
%                                     /          \
%                                    /         -------- |b>     State 2
%                  State 1    |a> -------
% 
% 0. syms wa wb we Wae wae Wbe wbe dae dbe real;
% 1. h=Hamiltonian(3);
% 2. h.addEnergies([wa,wb,we]);
% 3. h.addCoupling(1,3,Wae,wae);
% 4. h.addCoupling(2,3,Wbe,wbe);
% 5. h.defineStateDetuning(1,3,dae);
% 6. h.defineStateDetuning(2,3,dbe);
% 7. h.defineZero(we);
%
% Method I
%
% 8. h.unitaryTransformation();
%
% Method II
%
% 8. h.unitaryTransformation('Method','equations');
% 
% 
% 
%                               [0 0 0]
% 1. Creates a 3x3 empty matrix [0 0 0] and saves it to "hamiltonian" property. Creates also 3x3x2 empty matrix and saves it to
%                               [0 0 0]
% to the "coupling" property. Finally, creates empty vector [0, 0, 0] and saves it to "energies".
%                            [wa 0 0]
% 2. The hamiltonian becomes [0 wb 0], "energies" is now a vector [wa,wb,we].
%                            [0 0 we]                                                             [wa                  0   Wae/2*exp(i*wae*t)]
% 3. Adds a coupling between states 1 and 3. Hamiltonian obtains offdiagonal elements and becomes [0                   wb         0          ]. 
%                              [0 0 [wae Wae]]                                                    [Wae/2*exp(-i*wae*t) 0          we         ]
% Property "couplings" becomes [0 0     0    ]
%                              [0 0     0    ]                 [wa                  0                   Wae/2*exp(i*wae*t) ]
% 4. Different coupling is added. Hamiltonian becomes equal to [0                   wb                  Wbe/2*exp(-i*wbe*t)].
%                              [0 0 [wae Wae]]                 [Wae/2*exp(-i*wae*t) Wbe/2*exp(-i*wbe*t)          we        ]                                                                
% Property "couplings" becomes [0 0 [wbe Wbe]]
%                              [0 0     0    ] 
% 5. The detuning between states 1 and 3 is defined: wae=we-wa-dae. The function performs variable substitution.
% 6. The detuning between states 2 and 3 is defined: wbe=we-wb-dbe.
% 7. Sets excited state energy we to be zero and saves that information to 'zeroEnergy' property.
% 8. Unitary transformation is performed. In the Bloch equations, in order to perserve the form of the equations after 
% the transformation, the hamiltonian has to be transformed using the formula: Hnew=T'HT - iT'dT/dt, where A' is hermitian 
% conjugate. In general, the unitary transforamtion is performed by a matrix of the following form:
% 
%    [exp(i*a1*t)       ...0]
% T= [          \           ]
%    [0...       exp(i*an*t)]
% 
% After the transformation Hnew(j,k)=H(j,k)*exp(i(ak-aj)*t)+delta(j,k)*aj. We can then eliminate the time dependence by appropriate
% choice of exponents ai. We obtain a set of linear equations of the form w_jk=aj-ak, if the number of couplings (so non-zero off-diagonal elements)
% is at most equal to the number of states, which can be solved for all ai (Method II). In our example
%                                   [exp(i*wae*t)     0         0]
% the transformation matrix becomes [     0      exp(i*wae*t)   0] and is saved to property "transMatrix". The hamiltonian becomes
%                                   [     0           0         1]
% 
%        [we-dae   0     Wae/2]
% H_new= [ 0     we-dbe  Wbe/2] and is saved into property "transformed". Because we also specified the zero energy, substitution we=0 is performed  in the end.
%        [Wae/2  Wbe/2     we ]
%
% Method I uses a more general method. Here, all the diagonal elements of the unitary matrix are found by looking at 'coupling path' between the
% (i,i) state and the state defined as having zero energy (or a state close to it), if it's specified (otherwise state 1 energy is used). In our
% example, this method will lead to transformation matrix identical to what we obtained using Method II.
%



classdef Hamiltonian < handle
   properties
       hamiltonian
       transformed
       transMatrix
       energies
       couplings
       detunings
       stateGraph
       zeroEnergy
       brightStates
       darkStates
   end
   properties(Access=private)
       freqs
       cpl
       cplGr       
   end
   methods
       
       function obj=Hamiltonian(n)
           if isnumeric(n)
               n=round(n);
               H=sym('H',n);
               H(:,:)=sym(0);
               C=sym('C',[n n 2]);
               C(:,:,:,:)=sym(0);
               obj.hamiltonian=H;
               obj.energies=diag(H);
               obj.couplings=C;
               obj.detunings=[];
               obj.cpl=zeros(n,1);
               obj.freqs=[];
           else
               error('Dimension must be a numeric value (non-integers are rounded)')
           end
                
       end
       
       
       function obj=defineZero(obj,w)
          ens=obj.energies;
          if isempty(ens)
              error('You have to provide energies first')
          end
          
          envars=symvar(ens);
          
          if ~ismember(envars,w)
              error('Specified energy is not related to any of the energy levels')
          end
          
          obj.zeroEnergy=w; 
          
          obj.transformed=subs(obj.transformed,w,0);
           
       end
       
       
       function obj=addEnergies(obj,v)
           if length(v)==length(obj.hamiltonian)
               H=obj.hamiltonian;
               for i=1:length(v)
                   if isnumeric(v(i))
                       H(i,i)=sym(v(i));
                   else
                       H(i,i)=v(i);
                   end
               end
               obj.hamiltonian=H;
               obj.energies=diag(H);
           else
               error('Energy vector length must have the same number of elements as the diagonal of the hamiltonian (you must provide energies for all the states).')
           end
       end
       
       
       function obj=addCoupling(obj,S_i,S_f,W_r,w)
           H=obj.hamiltonian;
           C=obj.couplings;
           cp=obj.cpl;
           fr=obj.freqs;
         
           if isnumeric(S_i) && isnumeric(S_f)
               if S_i>length(H) || S_f>length(H)
                   error('Indecies cannot exceed size of the hamiltonian')
               end
               if S_i==S_f
                   error('Cannot couple with light the state with itself')
               end
               if isnumeric(W_r)
                   W_r=sym(W_r);
               end
               if isnumeric(w)
                   w=sym(w);
               end
               
               if ~ismember(w,fr)
                   fr=[fr,w];
               end
               C(S_i,S_f,1)=w;
               C(S_i,S_f,2)=W_r;


               syms t real; 
               H(S_i,S_f)=-W_r/2*exp(1i*w*t);
               H(S_f,S_i)=-conj(W_r)/2*exp(-1i*w*t);
               cp(S_i)=cp(S_i)+1;
               obj.cpl=cp;
               obj.freqs=fr;

           else
               error('Indecies have to be numeric values')
           end
           obj.hamiltonian=H;
           obj.couplings=C;
       end
       
       
       function obj=addPolyCoupling(obj,S_i,S_f,W_r,w)
           H=obj.hamiltonian;
           C=obj.couplings;
           cp=obj.cpl;
           fr=obj.freqs;
         
           if isnumeric(S_i) && isnumeric(S_f)
               if S_i>length(H) || S_f>length(H)
                   error('Indecies cannot exceed size of the hamiltonian')
               end
               if S_i==S_f
                   error('Cannot couple with light the state with itself')
               end
               if isnumeric(W_r)
                   W_r=sym(W_r);
               end
               if isnumeric(w)
                   w=sym(w);
               end
               
               if ~ismember(w,fr)
                   fr=[fr,w];
               end
               
               
%                if C(S_i,S_f,1,1)==0
                   C(S_i,S_f,1)=w;
                   C(S_i,S_f,2)=W_r;
%                else
%                    C(S_i,S_f,1,:)=[C(S_i,S_f,1,:end-1),w];
%                    C(S_i,S_f,2,:)=[C(S_i,S_f,2,:end-1),W_r];
%                end


               syms t real; 
               H(S_i,S_f)=H(S_i,S_f)-W_r/2*exp(1i*w*t);
               if imag(W_r)~=0
                   H(S_f,S_i)=H(S_f,S_i)-conj(W_r)/2*exp(-1i*w*t);
               else
                   H(S_f,S_i)=H(S_f,S_i)-W_r/2*exp(-1i*w*t);
               end
               cp(S_i)=cp(S_i)+1;
               obj.cpl=cp;
               obj.freqs=fr;

           else
               error('Indecies have to be numeric values')
           end
           obj.hamiltonian=H;
           obj.couplings=C;
       end
       
       
       function obj=createGraph(obj,L)
           import bloch.Dissipator

            if ~isa(L,'Dissipator')
                error('You have to provide the dissipator as an object of its own class')
            else
               C=obj.couplings;
               n=length(C);
               
               %Edges
               
               nodes_i=[];
               nodes_f=[];
               weights=[];
               labels={};
               
               for i=1:n
                   for j=1:n
                       if C(i,j,1)~=0
                           nodes_i=[nodes_i,i];
                           nodes_f=[nodes_f,j];
                           weights=[weights,1];
                           labels{end+1}=char(C(i,j,2));
                       end
                   end
               end

               
               Br=L.branching;
               Dr=L.decayR;
               vars=symvar(Dr);
               on=ones(1,length(vars));
               
               for i=1:n
                   for j=1:n
                       if Br(i,j)~=0
                           nodes_i=[nodes_i,i];
                           nodes_f=[nodes_f,j];
                           
                           w=subs(Br(i,j),vars,on);                          

                           weights=[weights,double(w)];
                           
                           dec=Br(i,j)*Dr(i);
                           
                           labels{end+1}=char(dec);
                       end
                   end
               end
               
               EdgeTable=table([nodes_i' nodes_f'],weights',labels','VariableNames',{'EndNodes', 'Weight', 'Label'});
               
               %Nodes
               
               names={};
               En=obj.energies;

               for j=1:length(En)
                    names{end+1}=char(En(j));
               end

               NodeTable=table(names','VariableNames',{'Label'});
               
               %Graph
               
               G=digraph(EdgeTable,NodeTable);

               obj.stateGraph=G; 
            end
       end

       
       function obj=plotGraph(obj)
           
           G=obj.stateGraph;
           figure
           plot(G,'NodeLabel',G.Nodes.Label,'EdgeLabel',G.Edges.Label,'Layout','force3')
           drawnow          
       end
       
       
       function obj=addSidebands(obj,center,order,depth,spacing)
           C=obj.couplings;
           H=obj.hamiltonian;
           Hf=obj.transformed;
           n=length(H);
           if length(Hf)~=n
               error('You can only add sidebands after you transform the hamiltonian');
           end
           
           t=sym('t','real');
           
           for ii=1:n
               for j=1:n
                   if C(ii,j,1)==center
                       A=H(ii,j);
                       B=H(j,ii);
                       Af=Hf(ii,j);
                       Bf=Hf(j,ii);
                       H(ii,j)=A*besselj(0,depth);
                       H(j,ii)=B*besselj(0,depth);
                       Hf(ii,j)=Af*besselj(0,depth);
                       Hf(j,ii)=Bf*besselj(0,depth);
                       for k=1:order
                           H(ii,j)=H(ii,j)-A*(besselj(k,depth)*exp(1i*k*spacing*t)+(-1)^k*besselj(k,depth)*exp(-1i*k*spacing*t));
                           H(j,ii)=H(j,ii)-B*(besselj(k,depth)*exp(-1i*k*spacing*t)+(-1)^k*besselj(k,depth)*exp(1i*k*spacing*t));
                           Hf(ii,j)=Hf(ii,j)-Af*(besselj(k,depth)*exp(1i*k*spacing*t)+(-1)^k*besselj(k,depth)*exp(-1i*k*spacing*t));
                           Hf(j,ii)=Hf(j,ii)-Bf*(besselj(k,depth)*exp(-1i*k*spacing*t)+(-1)^k*besselj(k,depth)*exp(1i*k*spacing*t));
                       end
                   end
               end
           end
           
           obj.hamiltonian=H;
           obj.transformed=Hf;
                          
       end
       
       
       function obj=addPhaseModulation(obj,center,depth,mod_freq)
           C=obj.couplings;
           H=obj.hamiltonian;
           Hf=obj.transformed;
           n=length(H);
           if length(Hf)~=n
               error('You can only add sidebands after you transform the hamiltonian');
           end
           
           t=sym('t','real');
           
           for ii=1:n
               for j=1:n
                   if C(ii,j,1)==center
                       A=H(ii,j);
                       B=H(j,ii);
                       Af=Hf(ii,j);
                       Bf=Hf(j,ii);
                       H(ii,j)=A*exp(1i*depth*sin(mod_freq*t));
                       H(j,ii)=B*exp(-1i*depth*sin(mod_freq*t));
                       Hf(ii,j)=Af*exp(1i*depth*sin(mod_freq*t));
                       Hf(j,ii)=Bf*exp(-1i*depth*sin(mod_freq*t));   
                   end
               end
           end
           
           obj.hamiltonian=H;
           obj.transformed=Hf;
                          
       end
       
       
       function obj=defineStateDetuning(obj,S_i,S_f,d)
           if isnumeric(d)
              d=sym(d);
           end

           H=obj.hamiltonian;
           Hf=obj.transformed;
           C=obj.couplings;
           D=obj.detunings;
           w0=obj.zeroEnergy;
           
           if ismember(d,D)
              error('Detuning has already been defined'); 
           end
           

           
           if C(S_i,S_f,1)~=0
               w=C(S_i,S_f,1);
           elseif C(S_f,S_i,1)~=0
               w=C(S_f,S_i,1);
           else
               error('There is no coupling between given states.')
           end
           
           
           Hf=subs(Hf,w,H(S_f,S_f)-H(S_i,S_i)-d);
           
           if ~isempty(w0)
            Hf=subs(Hf,w0,0);
           end
           
           obj.transformed=Hf;
           
         
           obj.detunings=[D;H(S_i,S_i),H(S_f,S_f),d,w];
        
       end
       
       
       function obj=defineEnergyDetuning(obj,E_i,E_f,d,w)
           if isnumeric(d)
              d=sym(d);
           end
    
           Hf=obj.transformed;
           D=obj.detunings;
           W=obj.freqs;
           w0=obj.zeroEnergy;
           
           if ismember(d,D)
              error('Detuning has already been defined'); 
           end
           
           if ~ismember(w,W)
               if ~ismember(w,symvar(W))
                error('Coupling with such frequency does not exist'); 
               end
           end

           

           Hf=subs(Hf,w,E_f-E_i-d);
           
           if ~isempty(w0)
            Hf=subs(Hf,w0,0);
           end
           
           obj.transformed=Hf;
           
           obj.detunings=[D;E_i,E_f,d,w];
           
       end
       
       
       function obj=unitaryTransformation(obj,varargin)
           
           defaultMethod='general';
           expectedMethods={'general','equations'};
           
           parser=inputParser;
           
           addParameter(parser,'Method',defaultMethod,@(x) any(validatestring(x,expectedMethods)));
           parse(parser,varargin{:})
           
           switch parser.Results.Method
               
               case 'equations'
                   obj.eqnTransform();
                   
               case 'general'
                   obj.generalTransform();
           end
       end
      
       
       function obj=changeBasis(obj,varargin)
           H=obj.transformed;
           ne=length(H);
           n=nargin-1;

           if n~=ne
               error('Number of vectors provided has to be equal to the number of states')
           end
           
           U=[];
           
           for i=1:n
               vec=varargin{i};

               if length(vec)~=ne
                
                   error('All vectors have to have the same length as number of states')
                   
               end
               s=size(vec);
               
               if (s(1)==1 && s(2)==1 && ne>1) || length(s)>2 || (s(1)>1 && s(2)>1)
                   error('You have to provide vectors.')
               end
               if s(1)==1
                   vec=vec.';
               end
               
               U=[U,vec];
               
           end
          
           obj.transformed=simplify(U'*H*U,'Steps',100);
                       
       end
       
       
       function obj=findDarkStates(obj,ground_states,excited_states)
           if ~isempty(intersect(ground_states,excited_states))
               disp('Indeces for ground and excited states have to be different')
               return
           end
           M=obj.transformed(ground_states,excited_states)';
           N=[];
           for i=1:length(excited_states)
               N=[N;M(i,:)/norm(squeeze(M(i,:)))];
           end
           N=simplify(N,'Steps',100);
           obj.brightStates=N;
           obj.darkStates=null(N);
           disp("There are "+num2str(length(ground_states)-rank(vpa(N)))+" dark states")
       end
       
       
       function Heff=adiabaticElimination(obj,ground_states,excited_states)
           if ~isempty(intersect(ground_states,excited_states))
               disp('Indeces for ground and excited states have to be different')
               return
           end
           n=length(obj.transformed);
           remainder=1:n;
           remainder=setdiff(remainder,[ground_states,excited_states]);
           P=zeros(n);
           Q=zeros(n);
           for i=ground_states
            P(i,i)=1;
           end
           for i=excited_states
            Q(i,i)=1;
           end
           D=diag(diag(obj.transformed));
           V=sym(zeros(n));
           V(ground_states,excited_states)=obj.transformed(ground_states,excited_states);
           V(excited_states,ground_states)=obj.transformed(excited_states,ground_states);
           invR=-(Q*D*Q+Q*V*Q);
%            disp(V)
%            disp(invR)
           R=V+V*(Q/invR)*V;
           
           Heff=P*D*P+P*R*P;
           
           if ~isempty(remainder)
              Heff(:,remainder)=obj.transformed(:,remainder);
              Heff(remainder,:)=obj.transformed(remainder,:);
               
           end
       end
   end
   
   
   methods(Access=private)
       
       function obj=eqnTransform(obj)
           
           w0=obj.zeroEnergy;
           H=obj.hamiltonian;
           C=obj.couplings;
           n=length(H);
           Eqns=[];
           t=sym('t','real');
           A=sym('a',[n 1]);
           
           try
               for i=1:n
                   for j=1:n
                        if C(i,j,1)~=0 && length(Eqns)<n
                            Eqns=[Eqns;C(i,j,1)==A(i)-A(j)]; %These are the linear equations that we need to solve to obtain transformation matrix T. 
                        end
                   end
               end
               sol=solve(Eqns,A);                         %Equations are solved.
           catch Er
               disp('Couldnt find unitary transformation. Using general method instead.')
               obj.generalTransform();
               return
           end
           
           T=sym('T',n);
           T(:,:)=sym(0);
           try
               for i=1:n
                   T(i,i)=exp(1i*sol.(['a' num2str(i)])*t);  %Solutions are puut into matrix T in the form exp(i*ai*t).
               end
           catch Er
               disp('Couldnt find unitary transformation. Using general method instead.')
               obj.generalTransform();
               return
           end
           
           
           H_f=T'*H*T-1i*T'*diff(T,t);
           H_f=simplify(H_f,'Steps',10);
           
          
           for i=1:n
               for j=i+1:n
                   if H_f(i,j)~=0
                       val=H_f(i,j)*2/C(i,j,2);  
                       if val~=1 && val~=-1
                           disp('Couldnt find unitary transformation. Using general method instead.')
                           obj.generalTransform();
                           return
                       end
                   end
               end
           end
           
           if ~isempty(obj.detunings)
               D=obj.detunings;
               for i=1:length(D(:,1))
                   w=D(i,4);
                   det=D(i,2)-D(i,1)-D(i,3);
                   H_f=subs(H_f,w,det);
               end
           end
           
           if ~isempty(w0)
            H_f=subs(H_f,w0,0);
           end
           
           obj.transMatrix=T;
           obj.transformed=H_f;
       end
    
       
       function obj=generalTransform(obj)

           H=obj.hamiltonian;
           if ~isempty(obj.transformed)
               H=obj.transformed;
           end
           w0=obj.zeroEnergy;
           if isempty(w0)
              disp('To use this method, you should specify the zero energy level. Use "defineZero" function. Now, the function will use H(1,1) as its baseline instead.')
              w0=H(1,1);
           end
           
           C=obj.couplings;
           cp=obj.cpl;
           n=length(H);
           t=sym('t','real');
           
           obj.couplingGraph();

           T=sym('T',n);
           T(:,:)=sym(0);
           for i=1:n
               T(i,i)=1;
           end
           for i=1:n
               phase=obj.shortestCouplingPath(i,w0);
               T(i,i)=T(i,i)*exp(1i*phase*t);                      
           end
           
           T=simplify(T);
           
           H_f=T'*H*T-1i*T'*diff(T,t);
           H_f=simplify(H_f,'Steps',10);
          


           if ~isempty(obj.detunings)
               D=obj.detunings;

               for i=1:length(D(:,1))
                   w=D(i,4);
                   det=(D(i,2)-D(i,1)-D(i,3));
                   H_f=subs(H_f,w,det);
               end
           end
           
           if ~isempty(w0)
            H_f=subs(H_f,w0,0);
           end
            
           obj.transMatrix=T;
           obj.transformed=H_f;
       end
       
       
       function obj=couplingGraph(obj)
               C=obj.couplings;
               n=length(C);
               
               %Edges
               
               nodes_i=[];
               nodes_f=[];
               weights=[];
               
               for i=1:n
                   for j=1:n
                       if C(i,j,1)~=0
                           nodes_i=[nodes_i,i];
                           nodes_f=[nodes_f,j];
                       end
                   end
               end
               
               for i=1:n
                   if ~ismember(i,nodes_i) && ~ismember(i,nodes_f)
                       nodes_i=[nodes_i,i];
                       nodes_f=[nodes_f,i];
                   end
               end
          
               
               EdgeTable=table([nodes_i' nodes_f'],'VariableNames',{'EndNodes'});
    
               %Graph
               
               G=graph(EdgeTable);


               
               obj.cplGr=G; 
            
       end
           
    
       function phase=shortestCouplingPath(obj,ind_i,w0)
           Gr=obj.cplGr;
           Ens=obj.energies;
           ind_w=0;   
           
           for i=1:length(Ens)
               s_en=symvar(Ens(i));

               if ismember(w0,s_en) || w0==Ens(i)
                   if ind_w>0
                    ind_w=[ind_w,i];
                   else
                       ind_w=i;
                   end
                  
               end

           end
           if ~ind_w
                 phase=0;
                 return
           end
           

           
           for ind=ind_w
%                if ismember(ind_i,ind_w)
%                    p=0;
%                    break
%                end
               p=shortestpath(Gr,ind_i,ind);

               if ~isempty(p)
                   break
               end
           end

           C=obj.couplings;
           phase=0;
           
           n=length(p);
           
           for j=1:n-1
               if C(p(j),p(j+1),1)~=0
                   phase=phase+C(p(j),p(j+1),1);
               else
                   phase=phase-C(p(j+1),p(j),1);
               end
           end
                
       end
       
   end
        
end


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
       




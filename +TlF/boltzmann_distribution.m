function[init]=boltzmann_distribution(rot_constant,temperature,J_g,varargin)   %Assumes 1Sigma+ state (L=0, S=0 -> Ja=0, Omega=0)
         
      defaultMethod='full';
      expectedMethods={'full','toy'};
           
      parser=inputParser;
           
      addParameter(parser,'Method',defaultMethod,@(x) any(validatestring(x,expectedMethods)));
      parse(parser,varargin{:})
           
      switch parser.Results.Method
               
        case 'toy'
              mult=2;
        case 'full'
              mult=4;
      end


    hbar=1.054*10^(-34); %[Js]
    k_b=1.381*10^(-23); % [J/K]
    Z=0;
    g=1;
    init=zeros(length(J_g),1);
    for i=1:length(J_g)
        g=mult*(2*J_g(i)+1);
        Z=Z+g*exp(-rot_constant*J_g(i)*(J_g(i)+1)*2*pi*hbar/(temperature*k_b));
    end
    

    for i=1:length(J_g)
        g=mult*(2*J_g(i)+1);
        init(i)=g*exp(-rot_constant*J_g(i)*(J_g(i)+1)*2*pi*hbar/(temperature*k_b))/Z;      
    end
end
function[G,E]=generateStates() 
    G=[];
    E=[];
    for n=1
        for s=[-1/2,1/2]
            j=n+s;
            for i=[-1/2,1/2]
                f=j+i;
                for m=-f:1:f
                    G=[G;[0,n,1/2,j,f,m]]; %[G/E,N,S,J,F,Mf]
                end
            end
        end
    end
    for j=1/2
        for i=[-1/2,1/2]
            f=j+i;
            for m=-f:1:f
                E=[E;[1,1,1/2,1/2,1/2,j,f,m,1]];%[G/E,lambda,S,sigma,omega,J,F,Mf,parity] (A^2Pi_1/2 state -> L=1, S=1/2, sigma=+-1/2, lambda=+-1, |omega|=1/2)
            end
        end
    end
end
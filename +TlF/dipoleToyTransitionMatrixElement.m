function[x]=dipoleToyTransitionMatrixElement(E,G,p)

    import TlF.*



    E_st=[]; 
    G_st=[]; 

            
              Je=E(2);
              We=E(3);
              Fe=E(4);
              Me=E(5);
              Pe=E(6);

              Jg=G(2);
              Wg=G(3);
              Fg=G(4);
              Mg=G(5);
              Pg=G(6);
              

              power=Fe-Me+Fg+2*Je-We+3/2;
              parg=(-1)^(Jg+Pg);
              pare=(-1)^(Je+Pe);
              del=1-delta_kr(parg,pare);
              sum=0;
              for q=[-1,0,1]
                  if We~=0 && Wg~=0
                    sum=sum+ThreeJSymbol(Je,-We,1,q,Jg,Wg)+(-1)^Pe*ThreeJSymbol(Je,We,1,q,Jg,Wg)+(-1)^Pg*ThreeJSymbol(Je,-We,1,q,Jg,-Wg)+(-1)^(Pg+Pe)*ThreeJSymbol(Je,We,1,q,Jg,-Wg);
                  elseif Wg==0 && We~=0
                    sum=sum+sqrt(2)*(ThreeJSymbol(Je,-We,1,q,Jg,Wg)+(-1)^Pe*ThreeJSymbol(Je,We,1,q,Jg,Wg));
                  elseif We==0 && Wg~=0
                    sum=sum+sqrt(2)*(ThreeJSymbol(Je,-We,1,q,Jg,Wg)+(-1)^Pg*ThreeJSymbol(Je,-We,1,q,Jg,-Wg));
                  else
                    sum=sum+2*ThreeJSymbol(Je,-We,1,q,Jg,Wg);
                  end  
              end
              x=((-1)^power)*del*sum*sqrt((2*Fg+1)*(2*Jg+1)*(2*Fe+1)*(2*Je+1))*ThreeJSymbol(Fe,-Me,1,p,Fg,Mg)*SixJSymbol(Jg,Fg,1/2,Fe,Je,1)*1/2;%/sqrt(factorial(Jg+3)*factorial(2-Jg));
           
    
end
    
    
 


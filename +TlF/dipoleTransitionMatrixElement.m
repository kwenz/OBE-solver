function[x]=dipoleTransitionMatrixElement(E,G,p,varargin)

    import TlF.*



    E_st=[]; 
    G_st=[]; 

    defaultMix='none';
    expectedMix={'none','include'};
           
      parser=inputParser;
           
      addParameter(parser,'Mixing',defaultMix,@(x) any(validatestring(x,expectedMix)));
      parse(parser,varargin{:});
           
      switch parser.Results.Mixing
               
        case 'none'
            
              Je=E(2);
              We=E(3);
              F1e=E(4);
              Fe=E(5);
              Me=E(6);
              Pe=E(7);

              Jg=G(2);
              Wg=G(3);
              F1g=G(4);
              Fg=G(5);
              Mg=G(6);
              Pg=G(7);
              

              power=Fe-Me+Fg+F1g+F1e+2*Je-We+3;
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
              x=((-1)^power)*del*sum*sqrt((2*Fg+1)*(2*Jg+1)*(2*Fe+1)*(2*Je+1)*(2*F1g+1)*(2*F1e+1))*ThreeJSymbol(Fe,-Me,1,p,Fg,Mg)*SixJSymbol(Jg,F1g,1/2,F1e,Je,1)*SixJSymbol(F1g,Fg,1/2,Fe,F1e,1)*1/2;%/sqrt(factorial(Jg+3)*factorial(2-Jg));
            
        case 'include'
              x=0;
      end
    
    
    
    
    
    
    
    
end
    
    
    
%     
%     
%     
%     if E(1)>0 || G(1)>0
%         Mixing_M=zeros(4,8); %|J,F1,F>: |1,1/2,0>,|1,1/2,1>,|1,3/2,1>,|1,3/2,2>,|2,3/2,1>,|2,3/2,2>,|2,5/2,2>,|3,5/2,2>
%         Mixing_M(1,:)=[1,0,0,0,0,0,0,0]; %E. Norrgard et al.
%         Mixing_M(2,:)=[0,0.9996,0.0203,0,-0.018,0,0,0];
%         Mixing_M(3,:)=[0,0.0267,-0.8518,0,0.5232,0,0,0];
%         Mixing_M(4,:)=[0,0,0,0.8482,0,-0.5294,-0.0138,0.0064];
%     end
%     
%     e_ind=0;
%     g_ind=0;
%     
%     if E(4)==1/2
%         if E(5)==0
%             e_ind=1;
%         else
%             e_ind=2;
%         end
%     else
%         if E(5)==1
%             e_ind=3;
%         else
%             e_ind=4;
%         end
%     end
%     
%     if G(4)==1/2
%         if G(5)==0
%             g_ind=1;
%         else
%             g_ind=2;
%         end
%     else
%         if G(5)==1
%             g_ind=3;
%         else
%             g_ind=4;
%         end
%     end
%   
%     
%     if E(1)>0      
%         for i=1:8
%             if Mixing_M(e_ind,i)~=0
%                 if i==1
%                     E_st=[E_st;[1,E(3),1/2,0,E(6),i]];
%                 elseif i==2
%                      E_st=[E_st;[1,E(3),1/2,1,E(6),i]];
%                 elseif i==3
%                      E_st=[E_st;[1,E(3),3/2,1,E(6),i]];
%                 elseif i==4
%                      E_st=[E_st;[1,E(3),3/2,2,E(6),i]];
%                 elseif i==5
%                      E_st=[E_st;[2,E(3),3/2,1,E(6),i]];
%                 elseif i==6
%                      E_st=[E_st;[2,E(3),3/2,2,E(6),i]];
%                 elseif i==7
%                      E_st=[E_st;[2,E(3),5/2,2,E(6),i]];
%                 else
%                      E_st=[E_st;[3,E(3),5/2,2,E(6),i]];                   
%                 end                    
%             end
%         end
%     else
%         E_st=[E(2),E(3),E(4),E(5),E(6)];
%     end
%     
%     if G(1)>0      
%         for i=1:8
%             if Mixing_M(g_ind,i)~=0
%                 if i==1
%                     G_st=[G_st;[1,G(3),1/2,0,G(6),i]];
%                 elseif i==2
%                      G_st=[G_st;[1,G(3),1/2,1,G(6),i]];
%                 elseif i==3
%                      G_st=[G_st;[1,G(3),3/2,1,G(6),i]];
%                 elseif i==4
%                      G_st=[G_st;[1,G(3),3/2,2,G(6),i]];
%                 elseif i==5
%                      G_st=[G_st;[2,G(3),3/2,1,G(6),i]];
%                 elseif i==6
%                      G_st=[G_st;[2,G(3),3/2,2,G(6),i]];
%                 elseif i==7
%                      G_st=[G_st;[2,G(3),5/2,2,G(6),i]];
%                 else
%                      G_st=[G_st;[3,G(3),5/2,2,G(6),i]];                   
%                 end                    
%             end
%         end 
%    
%     else
%         G_st=[G(2),G(3),G(4),G(5),G(6)];
%     end
% 
% 
% 
%     x=0;
%     for i=1:length(E_st(:,1))
%         for j=1:length(G_st(:,1))
%             Je=E_st(i,1);
%             We=E_st(i,2);
%             F1e=E_st(i,3);
%             Fe=E_st(i,4);
%             Me=E_st(i,5);
%             Pe=E(7);
% 
%             Jg=G_st(j,1);
%             Wg=G_st(j,2);
%             F1g=G_st(j,3);
%             Fg=G_st(j,4);
%             Mg=G_st(j,5);
%             Pg=G(7);
%             
%             partial_x=0;
%             
%             power=Fe-Me+Fg+F1g+F1e+2*Je-We+3;
%             for q=[-1,0,1]
%                 partial_x=partial_x+((-1)^power)*(1-delta(Pg,Pe))*sqrt((2*Fg+1)*(2*F1g+1)*(2*Jg+1)*(2*Fe+1)*(2*F1e+1)*(2*Je+1))*ThreeJSymbol(Fe,-Me,1,p,Fg,Mg)*ThreeJSymbol(Je,-We,1,q,Jg,Wg)*SixJSymbol(F1g,Fg,1/2,Fe,F1e,1)*SixJSymbol(Jg,F1g,1/2,F1e,Je,1);
%             end
%             
%             if E(1)>0 && G(1)==0
%                 x=x+partial_x*Mixing_M(e_ind,E_st(i,6));
%             elseif E(1)==0 && G(1)>0
%                 x=x+partial_x*Mixing_M(g_ind,G_st(j,6));
%             elseif E(1)>0 && G(1)>0
%                 x=x+partial_x*Mixing_M(g_ind,G_st(j,6))*Mixing_M(e_ind,E_st(i,6));
%             else
%                 x=partial_x;
%             end
%         end
%     end
% end



clear;
import bloch.*

warning('off','all');


%Hamiltonian

syms wa wb Dg we W1 W2 w1 w2 d1 d2 B ws real;


h=Hamiltonian(3);
h.addEnergies([wa,wb,we]);
h.addCoupling(1,3,W1,w1);
h.addCoupling(2,3,W2,w2);
h.defineStateDetuning(1,3,d1);
h.defineStateDetuning(2,3,d2);
h.defineZero(we);
h.unitaryTransformation();
h.addPhaseModulation(w1,B,ws);
% 
% 
disp(h.transformed);

h.changeBasis(1/sqrt(W1^2+W2^2)*[W1,W2,0],1/sqrt(W1^2+W2^2)*[W2,-W1,0],[0,0,1]);
disp(h.transformed);

return
% %Dissipator
% 
% syms G real;
% assume(G,'positive')
% L=Dissipator(3);
% L.addDecay(3,1,G/2);
% L.addDecay(3,2,G/2);
% 
% %Bloch equations
% 
% eq=BlochEqns(h,L);
% 
% 
% %Time evolution
% 
% % eq.necessaryVariables(); %Called first, to see which variables need to be substituted.
% % return
% % gamma=1;
% % % depth=5.3;
% % % spacing=1;
% % Om=10;
% % del=10/3;
% % Delta=10;
% % 
% % IC=zeros(3);  
% % IC(1,1)=1/2;
% % IC(2,2)=1/2;
% % 
% % eq.evolve(0,50,IC,[Delta,gamma,Om,del]);
% % eq.plotEvolution();
% % % 
% % pop=zeros(31,31);
% % 
% % parfor vard=1:31
% % 
% %     depth=0.3*(vard-1);
% %    
% %     for vars=1:31
% %         spacing=(vars-1)*0.3;
% %         eq.evolve(0,100,IC,[depth,Delta,gamma,Om,del,spacing]);
% %         pee=real(squeeze(eq.evolution(3,3,:)));
% %         lgt=round(4*length(pee(:,1))/5);
% %         res=mean(pee(lgt:end,1));
% %         pop(vard,vars)=res;
% %         
% %     
% %     end
% %     
% %     disp(vard)
% % end
% % 
% % surf(0:0.3:9,0:0.3:9,pop,'EdgeColor','none','facecolor','interp')
% % colormap('jet')
% % ylabel('Modulation depth')
% % xlabel('Spacing [\Gamma]')
% % title('Rabi rate = 1\Gamma');
% % c=colorbar;
% % c.Label.String='\rho_{ee}';
% % view(2)
% % drawnow        
% %         
%         
% 
% eq.solveSteady();
% % 
pee=eq.steadyState(3,3);
disp(pee)
% [num,den]=numden(pee);
% A=coeffs(den,d);
% B=coeffs(den,G);
% [C,D]=coeffs(den,[d,G]);
% disp(C)
% disp(D);




[res,res1]=getLorentzian(pee,d,G);
disp(res)
disp(res1)

% sol=solve(diff(pee,W)==0,W);
% 
% disp(sol)
% 
% pe=subs(pee,W,sol(2,1));
% pe=simplify(pe,'Steps',100);
% disp(pe)

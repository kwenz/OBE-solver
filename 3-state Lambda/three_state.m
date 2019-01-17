clear
import bloch.*

warning('off','all');


%Hamiltonian

syms w0 we D Wa Wb W wa w da db d w1 wb real


h=Hamiltonian(3);
h.addEnergies([w0,w0+D,we]);
h.addCoupling(1,3,W,w);
h.addCoupling(2,3,W,w);
% h.addCoupling(3,4,Wa,wa);
h.defineStateDetuning(1,3,d);
% h.defineStateDetuning(2,3,db);
% h.defineStateDetuning(3,4,da);
h.defineZero(we);
h.unitaryTransformation();
% h.addPhaseModulation(wae,B,wm)

% disp('Hamiltonian before and after unitary transformation')
% disp(h.hamiltonian)
% disp(h.transformed)
% h.changeBasis(1/sqrt(2)*[1,1,0],1/sqrt(2)*[1,-1,0],[0,0,1]);
% disp(h.transformed)

%%
%Dissipator
import bloch.*
syms G Ga real;
assume(G,'positive')
% assume(Ga,'positive')
L=Dissipator(3);
L.addDecay(3,1,G/2);
L.addDecay(3,2,G/2);
% U=1/sqrt(2)*[1,1,0;1,-1,0;0,0,sqrt(2)];
% L.addDecay(4,2,Ga/2);
% L.dissipator=simplify(U'*L.dissipator*U,'Steps',100);
% disp('Dissipator and branching ratios')
% disp(L.dissipator)
% disp(d.decayR)
% disp(d.branching)


IC=zeros(3);
IC(1,1)=1/2;
IC(2,2)=1/2;
%Bloch equations

eq=BlochEqns(h,L);
eq.necessaryVariables();
eq.evolve(0,50000,IC,[0.01,0.000000,0.0001,0]);
eq.plotEvolution()
    return
% %%
% eq.solveSteady();
% pab=eq.steadyState(1,1)-eq.steadyState(2,2);
% pab=simplify(subs(pab,[da,db,Wa,Wb,G],[-10*D,-9*D,W,W,1]),'Steps',200);
% solw=solve(diff(pab,W)==0,W);
% 
% %%
% disp(solw)
% %%
% figure
% 
%     fplot(solw(2))
%     xlim([0 0.5])
%     ylim([0 5])
%     hold on
%     syms x real;
%     fplot(x*4.555+0.06,'r --')
%     
% % % return
% % % 
% % % % dm=U'*dm*U;
% % %%
% % dm=subs(eq.steadyState(:,:),[da,db],[d,db-d]);
% % pdark=simplify((dm(1,1)+dm(2,2)-2*real(dm(1,2)))/2,'Steps',200);
% % pbright=simplify((dm(1,1)+dm(2,2)+2*real(dm(1,2)))/2,'Steps',200);
% % disp(limit(pdark,d,Inf))
% % disp(limit(dm(2,2),d,Inf))
% % return
% %%
% pd=subs(pdark,[G,d],[1,1]);
% pb=subs(pbright,[G,d],[1,1]);
% % disp(pd)
% % 
% % 
% figure
% fsurf(pd)
% ylim([0 5])
% 
% figure
% fsurf(pb)
% ylim([0 5])
% return
% % figure
% % fsurf((pb+pd)/2)
% % xlim([0 1])
% % return
% %%
% % Time evolution
% 
% eq.necessaryVariables(); %Called first, to see which variables need to be substituted.
% % return 
% 
% % 
% IC=zeros(3);  
% IC(1,1)=0;
% IC(2,2)=1;
% 
% gamma=1;
% 
% eq.evolve(0,50000,IC,[1,0.9,0.9,-10,-10.09]);
% eq.plotEvolution();
% return
% 
% 
% eq.eqnsRHS=subs(eq.eqnsRHS,[G,Wa,da],[1,0.1,0]);
% 
% for Delta=[0.02,2]
%    figure
%    k=1;
%    
%    for gamma_a=[0,0.0001,0.001,0.01,0.05,0.1,0.5,1,5]
%        
%         pop=zeros(51,51);
% 
%         parfor vard=1:51
% 
%             Omega=0.04*(vard-1);
%    
%             for vars=1:51
%                 det=(vars-26)*Delta/10;
% 
%                 eq.evolve(0,1000,IC,[Delta,gamma_a,Omega,det]);
%                 res=real(squeeze(eq.evolution(1,1,end)-eq.evolution(2,2,end)));
%                 pop(vard,vars)=res;
%             end
%             
%         end
%         subplot(3,3,k)
%         surf(-2.5*Delta:Delta/10:Delta*2.5,0:0.04:2,pop,'EdgeColor','none','facecolor','interp')
%         colormap('jet')
%         ylabel('\Omega [\Gamma]')
%         xlabel('\delta [\Gamma]')
%         c=colorbar;
%         c.Label.String='\rho_{11}-\rho_{22}';
%         view(2)
%         drawnow        
%         k=k+1;
%    end
%     
%     
% end
% % eq.evolve(0,50,IC,[Delta,gamma,Rabi,det]);
% % 
% % eq.plotEvolution();
% % return
% % 
% %
% %         
% 
% 
% % figure
% % plot(eq.evTime(1,:),squeeze(eq.evolution(3,3,:)))
% % xlabel('Time [1/\Gamma]')
% % ylabel('\rho_{ee}')
% 


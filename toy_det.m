clear;

import bloch.*

warning('off','all');


syms w0 D w1 w2 w3 we real;
syms W_L w_L W_m1 w_m1 W_m2 w_m2 W_m3 w_m3 real;
syms d_L d_m1 d_m2 d_m3 real;
syms Ga Gb real;
assume(Ga,'positive')
assume(Gb,'positive')





H=Hamiltonian(6);
H.addEnergies([w0,w0+D,w1,w2,w3,we]);
H.addCoupling(1,3,W_m1,w_m1);
H.addCoupling(2,3,W_m1,w_m1);
H.addCoupling(3,4,W_m2,w_m2);
H.addCoupling(4,5,W_m3,w_m3);
H.addCoupling(4,6,W_L,w_L);
H.defineZero(w0);
H.unitaryTransformation();
H.defineEnergyDetuning(w0,w1,d_m1,w_m1);
H.defineEnergyDetuning(w1,w2,d_m2,w_m2);
H.defineEnergyDetuning(w2,w3,d_m3,w_m3);
H.defineEnergyDetuning(w2,we,d_L,w_L);

disp(H.transformed)

L=Dissipator(6);
L.addDecay(6,4,Ga/2);
L.addDecay(6,2,Ga/4);
L.addDecay(6,1,Ga/4);
L.addDecay(3,2,Gb/4);
L.addDecay(3,1,Gb/4);

eq=BlochEqns(H,L);

IC=zeros(6);

IC(1,1)=1/10;
IC(2,2)=1/10;
IC(3,3)=1/5;
IC(4,4)=2/5;
IC(5,5)=1/5;

%%
gamma_a=1;
gamma_b=0;
% Omega_m1=0.01;
% Omega_m2=0.172;
% Omega_m3=0.0;
% Omega_L=1.46;
Delta=0.1;
% del_m1=1.1*Delta;
del_m2=0;
del_m3=0;
del_L=0;
% gambs=[0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,0.2,0.5,1,1.2,1.5,2,5,10];
% depth=2;
% modfr=3;
% 
% 
% figure
% k=0;
% for del_m1=[-Delta,0,Delta]
% % k=k+1;

%%
t_tot=1000;
no_switches=12;
t_step=floor(t_tot/no_switches);

%%

% 
% %
% 
% mat=[mat,real(eq.evolution(1,1,end))];
% % if k==1
% %     plot(mat(:,1),mat(:,2),'-ob')
% % elseif k==2
% %     plot(mat(:,1),mat(:,2),'-og')
% % else
% %     plot(mat(:,1),mat(:,2),'-or')
% % end
% % hold on
% % drawnow
% end
% 
% figure
% semilogx(gambs,mat,'b -o')
% 
% % end
% % 
% return
eq.eqnsRHS=subs(eq.eqnsRHS,[Ga,Gb,D,d_m2,d_m3,d_L],[gamma_a,gamma_b,Delta,0,0,0]);
%%
Params={W_L,0,1;W_m1,0,10*Delta;W_m2,0,1;W_m3,0,1;d_m1,0*Delta,10*Delta};
Sw_params=[1,0,1,1,1;1,1,0,0,1];

eq.optimizeParameters(0,t_tot,IC,Params,[1],'Integration','no','Popsize',15,'Iterations',10,'Popnumber',3,'Switching','yes','NoSwitches',no_switches,'SwitchingParameters',Sw_params);
%%

for sw=1:no_switches
    
    R_l=0.7;
    det=0.0942;
    
    if mod(sw,2)==1
        R_m1=0.0497;
        R_m2=0;
        R_m3=0;
    else
        R_m1=0;
        R_m2=0.53;
        R_m3=0.35;
    end
    

if sw==1
    eq.evolve(0,t_step,IC,[R_l,R_m1,R_m2,R_m3,det]);
else
    eq.intTime=[(sw-2)*t_step,(sw-1)*t_step];
    eq.lastSol=prevSol;
    eq.extendEvolution(t_step*sw,[R_l,R_m1,R_m2,R_m3,det]);
end

prevSol=eq.lastSol;
end

eq.plotEvolution()

return




figure
k=0;
for del_m1=[-2*Delta,-3*Delta,-4*Delta,-5*Delta,-6*Delta,-7*Delta]
pop=zeros(51,51);
%% 
k=k+1;
parfor vard=1:51

    Omega_m1=0.02*(vard-1);
   
    for vars=1:51
        Omega_m2=0.02*(vars-1);
%         Omega_m2=Omega_m1;
        eq.evolve(0,500,IC,[Delta,gamma_a,gamma_b,Omega_L,Omega_m1,Omega_m2,del_L,del_m1,del_m2]);          
        
%         pee=real(squeeze(eq.evolution(1,1,end)));
%         pee2=real(squeeze(eq.evolution(2,2,)));
%         lgt=round(4*length(pee(:,1))/5);
%         res=mean(pee(lgt:end,1))-mean(pee2(lgt:end,1));
        pop(vard,vars)=real(eq.evolution(1,1,end)-eq.evolution(2,2,end));
        
    
    end
    
    disp(vard)
end
%%
subplot(3,3,k)
surf(0:0.05:1,0:0.05:1,pop,'EdgeColor','none','facecolor','interp')
colormap('jet')
ylabel('\Omega_{\mu 1} [\Gamma]')
xlabel('\Omega_{\mu 2} [\Gamma]')
ylim([0 0.25])
c=colorbar;
caxis([-0.4 0.4])
c.Label.String='\rho_{11}-\rho_{22}';
view(2)
title(['\delta_{\mu 1} = ' num2str(del_m1/Delta) '\Delta'])
drawnow  
end
return

% % plot(eq.evTime(1,:),squeeze(eq.evolution(1,1,:)))
% % hold on
% % plot(eq.evTime(1,:),squeeze(eq.evolution(2,2,:)))
% % hold on
% % plot(eq.evTime(1,:),squeeze(eq.evolution(3,3,:)))
% % hold on
% % plot(eq.evTime(1,:),squeeze(eq.evolution(4,4,:)))
%%
eq.solveSteady();
% % 
% % 
% % p00=eq.steadyState(1,1);
% % pgg=eq.steadyState(2,2);
% % pff=eq.steadyState(3,3);
pee=eq.steadyState(1,1);
disp(pee)
return
% % 
% % res=subs(pee,[W_m,W_L,d_L,d_m,D,g],[Omega_m,Omega_L,del_L,del_m,Delta,gamma]);
% % disp(vpa(res))
% % 
% [re,re1]=getLorentzian(pee,d_L,g);
% disp(simplify(re))
% disp(re1)
%%
pe=subs(pee,[g,D,d_m,d_L],[gamma,0,del_m,del_L]);

% sol=solve(diff(pe,d_m)==0,d_m);
% disp(sol)

% % l=subs(sol,[W_L,W_m,d_m],[D,0.001,0]);
% disp(sol)
% fplot(sol,[-5 5])


% disp(subs(p00,[D,g,d_m,d_L,W_m],[1,1,0,0,0]));
% pe=simplify(subs(pe,[g,D],[1,4]),'Steps',100);
% pe=subs(pe,g,1);

disp(pe)

figure
h1=fsurf(pe, [0 20 0 20]);

h1.EdgeColor='none';
h1.FaceColor='interp';
colormap('jet')
c=colorbar;
c.Label.String='\rho_{ee}';
xlabel('\Omega_L [\Gamma]')
ylabel('\Omega_{\mu} [\Gamma]')
title('\Delta= 0')

view(2)

return
 pff=subs(pe,D,10);
 
 solaa=solve(diff(pff,a)==0,a);
 disp(solaa)
 figure
 fplot(solaa,[-2.5 2.5])
 
 
 
 pfff=simplify(subs(pe,b,solaa));
 disp(pfff)

 disp(limit(pfff,a,1/2))
solbb=solve(diff(pfff,a)==0);

figure
fsurf(pfff, [0 20 -1 0])
%   pff=limit(pfff,D,Inf);
 
re=subs(re,[W_L,d_m],[D,-D/4]);
disp(re)

return




return

ans=(4*D^6 - 8*D^5*d_L + 16*D^5*d_m + 2*D^4*W_L^2 + 2*D^4*W_m^2 + 12*D^4*d_L^2 - 24*D^4*d_L*d_m + 20*D^4*d_m^2 + 3*D^4*g^2 + 8*D^3*W_L^2*d_m - 4*D^3*W_m^2*d_L + 4*D^3*W_m^2*d_m + 40*D^3*d_L^2*d_m - 16*D^3*d_L*d_m^2 + 8*D^3*d_m^3 + 10*D^3*d_m*g^2 + 2*D^2*W_L^4 + 4*D^2*W_L^2*W_m^2 + 8*D^2*W_L^2*d_m^2 + 2*D^2*W_m^4 + 8*D^2*W_m^2*d_L^2 + 8*D^2*W_m^2*d_L*d_m + 8*D^2*W_m^2*d_m^2 + 2*D^2*W_m^2*g^2 + 40*D^2*d_L^2*d_m^2 + 16*D^2*d_L*d_m^3 + 8*D^2*d_m^4 + 10*D^2*d_m^2*g^2 + 8*D*W_L^4*d_m - 8*D*W_L^2*W_m^2*d_L + 8*D*W_L^2*W_m^2*d_m - 8*D*W_m^4*d_L + 8*W_L^4*d_m^2 - 16*W_L^2*W_m^2*d_L*d_m + 8*W_m^4*d_L^2 + 2*W_m^4*g^2);
an=subs(ans,W_m,0);
disp(simplify(ans-an,'Steps',100))



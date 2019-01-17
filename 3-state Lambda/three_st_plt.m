% % plot(solution.x(1,:),solution.y(1,:))
% % hold on
% % plot(solution.x(1,:),solution.y(5,:))
% % hold on
% % plot(solution.x(1,:),solution.y(9,:))
% % % xlim([0 2])
% % legend('State A', 'State B', 'State E')
% disp(simplify(sol.Y5+sol.Y9,'Steps',100));
% % 
% % pope1=subs(sol.Y5,[G,D,d],[1,0.2,0.1]);
% % pope2=subs(sol.Y9,[G,D,d],[1,0.2,0.1]);
% limit(sol.Y5+sol.Y9,W_L,Inf)
% % fplot(pope1+pope2,[0 4])
% % p_12=[];
% % for x=real(Doublet(1,:))
% %     p_12=[p_12,(1i*Rabi_ae/2)/(2*Gamma_a/2+1i*x-1i*Rabi_be^2/(4*x))];
% % end
% % disp(imag(Doublet(2,:)))
% % 
% % 
% % subplot(2,2,1)
% % plot(real(Doublet(1,:)/(2*Gamma_a)),imag(Doublet(2,:)*2*Gamma_a/Rabi_ae))
% % ylim([0 1]);
% % xlabel('Probe detuning [\Gamma]')
% % ylabel('Im(\Gamma \rho_{ae}/\omega_p)')
% % title('Numerical result in weak probe limit')
% % 
% % subplot(2,2,2)
% % plot(real(Doublet(1,:)/(2*Gamma_a)),imag(p_12*2*Gamma_a/Rabi_ae))
% % xlabel('Probe detuning [\Gamma]')
% % ylabel('Im(\Gamma \rho_{ae}/\omega_p)')
% % title('Analytical result in weak probe limit')
% % 
% % subplot(2,2,3)
% % plot(real(Doublet(1,:)/(2*Gamma_a)),real(Doublet(2,:)*2*Gamma_a/Rabi_ae))
% % xlabel('Probe detuning [\Gamma]')
% % ylabel('Re(\Gamma \rho_{ae}/\omega_p)')
% % title('Numerical result in weak probe limit')
% % 
% % subplot(2,2,4)
% % plot(real(Doublet(1,:)/(2*Gamma_a)),real(p_12*2*Gamma_a/Rabi_ae))
% % xlabel('Probe detuning [\Gamma]')
% % ylabel('Re(\Gamma \rho_{ae}/\omega_p)')
% % title('Analytical result in weak probe limit')
h=- D/12 - ((13*D^2)/144 + W_L^2/6 + W_m^2/6)/(2*(- (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (D*W_L^2)/16 + ((- (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (D*W_L^2)/16 + (53*D^3)/1728)^2 - ((13*D^2)/144 + W_L^2/6 + W_m^2/6)^3)^(1/2) + (53*D^3)/1728)^(1/3)) - ((D*W_L^2)/16 - (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (((D*W_L^2)/16 - (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (53*D^3)/1728)^2 - ((13*D^2)/144 + W_L^2/6 + W_m^2/6)^3)^(1/2) + (53*D^3)/1728)^(1/3)/2 + (3^(1/2)*(((13*D^2)/144 + W_L^2/6 + W_m^2/6)/(- (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (D*W_L^2)/16 + ((- (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (D*W_L^2)/16 + (53*D^3)/1728)^2 - ((13*D^2)/144 + W_L^2/6 + W_m^2/6)^3)^(1/2) + (53*D^3)/1728)^(1/3) - ((D*W_L^2)/16 - (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (((D*W_L^2)/16 - (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (53*D^3)/1728)^2 - ((13*D^2)/144 + W_L^2/6 + W_m^2/6)^3)^(1/2) + (53*D^3)/1728)^(1/3))*1i)/2;
g=- D/12 - ((13*D^2)/144 + W_L^2/6 + W_m^2/6)/(2*(- (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (D*W_L^2)/16 + ((- (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (D*W_L^2)/16 + (53*D^3)/1728)^2 - ((13*D^2)/144 + W_L^2/6 + W_m^2/6)^3)^(1/2) + (53*D^3)/1728)^(1/3)) - ((D*W_L^2)/16 - (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (((D*W_L^2)/16 - (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (53*D^3)/1728)^2 - ((13*D^2)/144 + W_L^2/6 + W_m^2/6)^3)^(1/2) + (53*D^3)/1728)^(1/3)/2 - (3^(1/2)*(((13*D^2)/144 + W_L^2/6 + W_m^2/6)/(- (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (D*W_L^2)/16 + ((- (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (D*W_L^2)/16 + (53*D^3)/1728)^2 - ((13*D^2)/144 + W_L^2/6 + W_m^2/6)^3)^(1/2) + (53*D^3)/1728)^(1/3) - ((D*W_L^2)/16 - (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (((D*W_L^2)/16 - (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (53*D^3)/1728)^2 - ((13*D^2)/144 + W_L^2/6 + W_m^2/6)^3)^(1/2) + (53*D^3)/1728)^(1/3))*1i)/2;
f=((13*D^2)/144 + W_L^2/6 + W_m^2/6)/(- (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (D*W_L^2)/16 + ((- (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (D*W_L^2)/16 + (53*D^3)/1728)^2 - ((13*D^2)/144 + W_L^2/6 + W_m^2/6)^3)^(1/2) + (53*D^3)/1728)^(1/3) - D/12 + ((D*W_L^2)/16 - (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (((D*W_L^2)/16 - (D*(D^2/4 + W_L^2/2 + W_m^2/2))/24 + (53*D^3)/1728)^2 - ((13*D^2)/144 + W_L^2/6 + W_m^2/6)^3)^(1/2) + (53*D^3)/1728)^(1/3);
f=subs(f,[D,W_L],[5,5]);
f=real(f);
g=subs(g,[D,W_L],[5,5]);
g=real(g);
h=subs(h,[D,W_L],[5,5]);
h=real(h);


figure
fplot(f,[0 50])
hold on
fplot(g,[0 50])
hold on
fplot(h,[0 50])


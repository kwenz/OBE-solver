function [gammap,lor]=lorentzianForm(poly,delta,gamma)
    syms Gp real;
    
    [num,den]=numden(poly);
    r=num;
    
    
    C=coeffs(den,[delta,gamma],'All');
 


    a=C(1,3);
    b=C(2,3);
    c=C(3,1);
    d=C(3,3);
    

%     
    gammap=simplify(2*sqrt(c/a*gamma^2+d/a-b^2/(4*a^2)),'Steps',100);
    lor=(r/a/((delta+b/(2*a))^2+Gp^2/4));
    
end
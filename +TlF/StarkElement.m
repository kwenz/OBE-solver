function[x]=StarkElement(E,G,p,varargin)

    import TlF.*
    defaultSt='coupled';
    expectedSt={'coupled','uncoupled'};
           
      parser=inputParser;
           
      addParameter(parser,'Scheme',defaultSt,@(x) any(validatestring(x,expectedSt)));
      parse(parser,varargin{:});
           
      switch parser.Results.Scheme
               
        case 'coupled'
            
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
              

              power=Fe-Me+Fg+F1g+F1e+2*Je+3;
              x=((-1)^power)*sqrt((2*Fg+1)*(2*Jg+1)*(2*Fe+1)*(2*Je+1)*(2*F1g+1)*(2*F1e+1))*ThreeJSymbol(Fe,-Me,1,p,Fg,Mg)*SixJSymbol(Jg,F1g,1/2,F1e,Je,1)*SixJSymbol(F1g,Fg,1/2,Fe,F1e,1)*ThreeJSymbol(Je,0,1,0,Jg,0);
       case 'uncoupled'
               
              Je=E(1);
              mje=E(2);
              m1e=E(4);
              m2e=E(6);

              Jg=G(1);
              mjg=G(2);
              m1g=G(4);
              m2g=G(6);
              if Je==Jg && mjg==mje && m1g==m1e && m2g==m2e && Je>0 && p==0
                  x=(Je*(Je+1)-3*mje^2)/(Je*(Je+1)*(2*Je-1)*(2*Je+3));
              else
                power=2*Je-mje;
                x=((-1)^power)*sqrt((2*Je+1)*(2*Jg+1))*delta_kr(m1g,m1e)*delta_kr(m2g,m2e)*ThreeJSymbol(Je,mje,1,p,Jg,mjg)*ThreeJSymbol(Je,0,1,0,Jg,0);
              end
      end
end
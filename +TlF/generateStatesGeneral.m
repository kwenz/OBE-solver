function[states]=generateStates(J_g,Wg,Pg,J_e,We,Pe) %P=0 - e, P=1 - f (parity type of the excited state); Ground -> X1Sigma+, Excited B3Pi1
    G=[];
    E=[];
    for j=J_g
       if j==0
          G=[[0,0,Wg,1/2,0,0,Pg];[0,0,Wg,1/2,1,-1,Pg];[0,0,Wg,1/2,1,0,Pg];[0,0,Wg,1/2,1,1,Pg]];  %(ground (0) vs excited (1),J,W,F1,F,mF,parity type) 
       else
            for F1=[-1/2,1/2]
               for F2=[-1/2,1/2]
                      F=j+F1+F2;
                      for mF=-F:1:F
                          G=[G;[0,j,Wg,j+F1,F,mF,Pg]];
                      end
                  end
              end
          end
       end
    end
    for j=J_e
       if j==0
          E=[[1,0,We,1/2,0,0,P];[1,0,We,1/2,1,-1,Pe];[1,0,We,1/2,1,0,Pe];[1,0,We,1/2,1,1,Pe]];  %(X or B, J,|W|,F1,F,mF,parity) 
       else

              for F1=[-1/2,1/2]
                  for F2=[-1/2,1/2]
                      F=j+F1+F2;
                      for mF=-F:1:F
                          E=[E;[1,j,We,j+F1,F,mF,Pe]];
                      end
                  end
              end

       end
    end
    states=[G;E];
end
function[states]=generateStates(J_g,J_e,P) %P=0 - e, P=1 - f (parity type of the excited state); Ground -> X1Sigma+, Excited B3Pi1
    G=[];
    E=[];
    for j=J_g
       if j==0
          G=[[0,0,0,1/2,0,0,0];[0,0,0,1/2,1,-1,0];[0,0,0,1/2,1,0,0];[0,0,0,1/2,1,1,0]];  %(X or B,J,W,F1,F,mF,parity type) 
       else
          for W=0
              for F1=[-1/2,1/2]
                  for F2=[-1/2,1/2]
                      F=j+F1+F2;
                      for mF=-F:1:F
                          G=[G;[0,j,W,j+F1,F,mF,0]];
                      end
                  end
              end
          end
       end
    end
    for j=J_e
       if j==0
          E=[[1,0,1,1/2,0,0,P];[1,0,1,1/2,1,-1,P];[1,0,1,1/2,1,0,P];[1,0,1,1/2,1,1,P]];  %(X or B, J,|W|,F1,F,mF,parity) 
       else
          for W=1
              for F1=[-1/2,1/2]
                  for F2=[-1/2,1/2]
                      F=j+F1+F2;
                      for mF=-F:1:F
                          E=[E;[1,j,W,j+F1,F,mF,P]];
                      end
                  end
              end
          end
       end
    end
    states=[G;E];
end
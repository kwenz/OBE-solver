function[states]=generateDecoupledStates(J_g) %only ground states
    G=[];
    for j=0:1:J_g     %(J,MJ,I1,MI1,I2,MI2) 
          for MJ=-j:1:j
             for MI1=[-1/2,1/2]
                  for MI2=[-1/2,1/2]
                 
                          G=[G;[j,MJ,1/2,MI1,1/2,MI2]];
                  end
              end
          end
       
    end

    states=G;
end
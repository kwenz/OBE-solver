function x=GramSchmidt(varargin)
    n=nargin;
    x=[];

    for i=1:n
        vec=varargin{i};
        if length(vec)~=n
             error('All vectors have to have the same length and be equal to the number of vectors given.')
        end
        s=size(vec);
               
        if (s(1)==1 && s(2)==1 && n>1) || length(s)>2 || (s(1)>1 && s(2)>1)
           error('You have to provide vectors.')
        end
        if s(1)==1
           vec=vec.';
        end
        
        if i==1
            vec=vec/norm(vec);
            for ind=1:length(vec)
                if isnumeric(vec(ind))
                if abs(vec(ind))<10^(-6)
                    vec(ind)=0;
                end
                end
            end
            x=[x,vec];
        else
            v_gm=vec;
            if dot(v_gm,x(:,1))==0 && i==2
                x=[x,v_gm];
                disp(v_gm)
                continue
            end
            
            for j=1:i-1
                xv=x(:,j);
                v_gm=v_gm-dot(xv,vec)/dot(xv,xv)*xv;
            end
            v_gm=v_gm/norm(v_gm);
            for ind=1:length(v_gm)
                if isnumeric(v_gm(ind))
                if abs(v_gm(ind))<10^(-6)
                    v_gm(ind)=0;
                end
                end
            end
            if ~isnumeric(v_gm)
                v_gm=simplify(v_gm,'Steps',200);
            end
            x=[x,v_gm];
            
%             disp('Done')
            disp(v_gm)
        end
               

    end
    

    
    if ~isnumeric(x)    
        x=simplify(x,'Steps',300);
    end

end
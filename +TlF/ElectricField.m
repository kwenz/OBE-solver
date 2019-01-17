classdef ElectricField < handle
    
    properties
        direction
        Einternal
        Ecartesian
        Ecircular
        coords
    end
 
    
    methods
        
        function obj=ElectricField(amplitudes,direction)
            if length(amplitudes)~=3 || length(direction)~=3
                error('Please provide full vectors with three components')
            end
            
            if amplitudes(3)~=0
                error('There can be no electric field in the propagation direction')
            end
            
            obj.direction=direction;
            obj.Einternal=amplitudes;
        end

        function obj=getCoords(obj)
            vec=obj.direction;
            k1=vec(1);
            k2=vec(2);
            k3=vec(3);
            
            if k2==0 && k3==0
                xp=[0,1,0];
                yp=[0,0,1];

            elseif k1==0 && k3==0
                xp=[0,0,1];
                yp=[1,0,0];
                
            elseif k1==0 && k2==0
                xp=[1,0,0];
                yp=[0,1,0];
            else

                xp=[-k2,k1,0];
                if k3~=0 
                    yp=[k1,k2,-(k1^2+k2^2)/k3];
                else
                    yp=[0,0,-1];
                end
            end
            
            zp=vec;
                        
            xp=xp/norm(xp);
            yp=yp/norm(yp);
            zp=zp/norm(zp);

            x=[xp;yp;zp];

            obj.coords=x;
        end


        function obj=transformEfield(obj)

            C=obj.coords;
            xp=C(1,:);
            yp=C(2,:);
            zp=C(3,:);
            
            vec=obj.Einternal;

            Ex=vec(1);
            Ey=vec(2);


            theta=acos(zp(3));
            if zp(1)==0 && zp(2)==0
                phi=0;
            elseif zp(1)==0 && zp(2)>0
                phi=-pi/2;
            elseif zp(1)==0 && zp(2)<0
                phi=pi/2;
            elseif zp(2)==0 && zp(1)>0
                phi=0;
            elseif zp(2)==0 && zp(1)<0
                phi=pi;
            else                
                phi=atan(zp(2)/-zp(1));
            end
            ct=cos(theta);
            st=sin(theta);
            if abs(ct)<10^(-6)
                if theta<pi
                    ct=0;
                    st=1;
                else
                    ct=0;
                    st=-1;
                end
                
            end
 
            ep=exp(1i*phi);
            if abs(real(ep))<10^(-6)
                ep=imag(ep);
            elseif abs(imag(ep))<10^(-6)
                ep=real(ep);
            end

            Ep=simplify((Ex+1i*Ey*ct)/sqrt(sym('2'))*ep);
            Em=simplify((Ex-1i*Ey*ct)/sqrt(sym('2'))*ep);
            Ez=simplify(-Ey*st);

            E=[Em,Ez,Ep];

            obj.Ecircular=E;
            
            Ecx=simplify((Ep+Em)/sqrt(sym('2')));
            Ecy=simplify(1i*(Ep-Em)/sqrt(sym('2')));
            

            obj.Ecartesian=[Ecx,Ecy,Ez];


        end
    end
end

% 
% 
% classdef ElectricField < handle
%     
%     properties
%         direction
%         Einternal
%         Ecartesian
%         Ecircular
%         coords
%     end
%  
%     
%     methods
%         
%         function obj=ElectricField(amplitudes,direction)
%             if length(amplitudes)~=3 || length(direction)~=3
%                 error('Please provide full vectors with three components')
%             end

%             
%             obj.direction=direction;
%             obj.Ecartesian=amplitudes;
%         end
% 
%         function obj=getCoords(obj)
%             vec=obj.direction;
%             k1=vec(1);
%             k2=vec(2);
%             k3=vec(3);
%             
%             if k2==0 && k3==0
%                 xp=[0,1,0];
%                 yp=[0,0,1];
% 
%             elseif k1==0 && k3==0
%                 xp=[0,0,1];
%                 yp=[1,0,0];
%                 
%             elseif k1==0 && k2==0
%                 xp=[1,0,0];
%                 yp=[0,1,0];
%             else
% 
%                 xp=[-k2,k1,0];
%                 if k3~=0 
%                     yp=[k1,k2,-(k1^2+k2^2)/k3];
%                 else
%                     yp=[0,0,-1];
%                 end
%             end
%             
%             zp=vec;
%                         
%             xp=xp/norm(xp);
%             yp=yp/norm(yp);
%             zp=zp/norm(zp);
% 
%             x=[xp;yp;zp];
% 
%             obj.coords=x;
%         end
% 
% 
%         function obj=transformEfield(obj)
% 
%             C=obj.coords;
%             xp=C(1,:);
%             yp=C(2,:);
%             zp=C(3,:);
%             
%             vec=obj.Ecartesian;
% 
%             Ex=simplify(dot(vec,xp));
%             Ey=simplify(dot(vec,yp));
% 
%             if simplify(dot(vec,zp))~=0
%                 error('Field not perpendicular to propagation direction')
%             end
% 
%             theta=acos(zp(3));
%             if zp(1)==0 && zp(2)==0
%                 phi=0;
%             elseif zp(1)==0 && zp(2)>0
%                 phi=-pi/2;
%             elseif zp(1)==0 && zp(2)<0
%                 phi=pi/2;
%             elseif zp(2)==0 && zp(1)>0
%                 phi=0;
%             elseif zp(2)==0 && zp(1)<0
%                 phi=pi;
%             else                
%                 phi=atan(zp(2)/-zp(1));
%             end
%             ct=cos(theta);
%             st=sin(theta);
%             if abs(ct)<10^(-6)
%                 if theta<pi
%                     ct=0;
%                     st=1;
%                 else
%                     ct=0;
%                     st=-1;
%                 end
%                 
%             end
%  
%             ep=exp(1i*phi);
%             if abs(real(ep))<10^(-6)
%                 ep=imag(ep);
%             elseif abs(imag(ep))<10^(-6)
%                 ep=real(ep);
%             end
% 
%             Ep=simplify((Ex+1i*Ey*ct)/sqrt(sym('2'))*ep);
%             Em=simplify((Ex-1i*Ey*ct)/sqrt(sym('2'))*ep);
%             Ez=simplify(-Ey*st);
% 
%             E=[Em,Ez,Ep];
% 
%             obj.Ecircular=E;
% 
%             obj.Einternal=[Ex,Ey,0];
% 
% 
%         end
%     end
% end
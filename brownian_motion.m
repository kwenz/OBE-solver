function res=brownian_motion(N,T,mean,st_dev)
    vec=1/sqrt(N-1/T).*(st_dev.*randn(N-1,1)+mean);
    for i=2:N-1
        vec(i)=vec(i)+vec(i-1);
    end
    DT=0:T/(N-1):T;
%     disp(size(vec));
%     disp(size(DT));
    vec=[0;vec];
    res=[vec,DT'];
end
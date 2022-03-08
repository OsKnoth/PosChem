function [TOut,YOut] = RungeKuttaPat3Method(tBegin,tEnd,y,dt,RK)
global Reak Species dtMax
t=tBegin;
% YOut=zeros(size(y,1),(tEnd-tBegin)/dt);
% TOut=zeros(1,(tEnd-tBegin)/dt);
TOut(1)=t;
YOut(:,1)=y;

Y=zeros(size(y,1),RK.stage);
MY(RK.stage).M=sparse(size(y,1),size(y,1));
% Sigma=zeros(size(y,1),1);
p=3*RK.a(2,1)*(RK.a(3,1)+RK.a(3,2))*RK.b(3);

i=1;
TOut(i)=t;
YOut(:,i)=y;
while t<tEnd
    Y(:,1)=y;%Y1
    MY(1).M=MHatMatrix(Reak,Species,t,Y(:,1));
    M=RK.a(2,1)*MY(1).M;
    for k=1:size(y,1)
        M(:,k)=M(:,k)/max(Y(k,1),eps);
    end
    Y(:,2)=(speye(size(y,1))-dt*M)\y;
    
    MY(2).M=MHatMatrix(Reak,Species,t,Y(:,2));
    M=RK.a(3,1)*MY(1).M+RK.a(3,2)*MY(2).M;
    for k=1:size(y,1)
        rho=y(k)*(Y(k,2)/y(k))^(1/p);%zu Y(3)
        M(:,k)=M(:,k)/max(rho,eps);
    end
    Y(:,3)=(speye(size(y,1))-dt*M)\y;
    
    M=RK.beta(1)*MY(1).M+RK.beta(2)*MY(2).M;
    for j=1:size(y,1)
        my=y(j)*(Y(j,2)/y(j))^(1/RK.a(2,1));
        M(:,j)=M(:,j)/max(my,eps);
    end
    Sigma=(speye(size(y,1))-dt*M)\y;
    
    MY(3).M=MHatMatrix(Reak,Species,t,Y(:,3));
    M=RK.b(1)*MY(1).M+RK.b(2)*MY(2).M+RK.b(3)*MY(3).M;
    for j=1:size(y,1)
        M(:,j)=M(:,j)/max(Sigma(j,1),eps);
    end
    y=(speye(size(y,1))-dt*M)\y;
    t=t+dt;
    i=i+1;
    TOut(i+1)=t;
    t
    YOut(:,i+1)=y;
    dt=min(2*dt,dtMax);
    dt=min(dt,tEnd-t);
end
TOut=TOut';
TOut(end)
YOut=YOut';
end


function [TOut,YOut] = RungeKuttaPat3PD(tBegin,tEnd,y,dt,RK)
global dtMax
t=tBegin;
eps=1.e-10;
% YOut=zeros(size(y,1),(tEnd-tBegin)/dt);
% TOut=zeros(1,(tEnd-tBegin)/dt);
TOut(1)=t;
YOut(:,1)=y;

Y=zeros(size(y,1),RK.stage);
PY(RK.stage).M=sparse(size(y,1),size(y,1));
DY(RK.stage).M=sparse(size(y,1),size(y,1));
% Sigma=zeros(size(y,1),1);
p=3*RK.a(2,1)*(RK.a(3,1)+RK.a(3,2))*RK.b(3);

i=1;
TOut(i)=t;
YOut(:,i)=y;
while t<tEnd
    Y(:,1)=y;%Y1
    [PY(1).M,DY(1).M]=PDMatrix(t,Y(:,1));
    P=RK.a(2,1)*PY(1).M;
    D=RK.a(2,1)*DY(1).M;
    for k=1:size(y,1)
        D(k,:)=D(k,:)/max(Y(k,1),eps);
        P(:,k)=P(:,k)/max(Y(k,1),eps);
    end
    Y(:,2)=(speye(size(y,1))-dt*(P-D))\y;
    
    [PY(2).M,DY(2).M]=PDMatrix(t,Y(:,2));
    P=RK.a(3,1)*PY(1).M+RK.a(3,2)*PY(2).M;
    D=RK.a(3,1)*DY(1).M+RK.a(3,2)*DY(2).M;
    for k=1:size(y,1)
        rho=y(k)*(Y(k,2)/y(k))^(1/p);%zu Y(3)
        D(k,:)=D(k,:)/max(rho,eps);
        P(:,k)=P(:,k)/max(rho,eps);
    end
    Y(:,3)=(speye(size(y,1))-dt*(P-D))\y;
    
    P=RK.beta(1)*PY(1).M+RK.beta(2)*PY(2).M;
    D=RK.beta(1)*DY(1).M+RK.beta(2)*DY(2).M;
    for k=1:size(y,1)
        my=y(k)*(Y(k,2)/y(k))^(1/RK.a(2,1));
        D(k,:)=D(k,:)/max(my,eps);
        P(:,k)=P(:,k)/max(my,eps);
    end
    Sigma=(speye(size(y,1))-dt*(P-D))\y;
    
    
    [PY(3).M,DY(3).M]=PDMatrix(t,Y(:,3));
    P=RK.b(1)*PY(1).M+RK.b(2)*PY(2).M+RK.b(3)*PY(3).M;
    D=RK.b(1)*DY(1).M+RK.b(2)*DY(2).M+RK.b(3)*DY(3).M;
    for k=1:size(y,1)
        D(k,:)=D(k,:)/max(Sigma(k,1),eps);
        P(:,k)=P(:,k)/max(Sigma(k,1),eps);
    end
    y=(speye(size(y,1))-dt*(P-D))\y;
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


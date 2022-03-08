function [TOut,YOut] = RungeKuttaExpMethod(tBegin,tEnd,y,dt,RK)
global Reak Species dtMax
t=tBegin;
% YOut=zeros(size(y,1),2*(tEnd-tBegin)/dtMax);
% TOut=zeros(1,2*(tEnd-tBegin)/dtMax);
TOut(1)=t;
YOut(:,1)=y;

Y=zeros(size(y,1),RK.stage);
MY(RK.stage).M=sparse(size(y,1),size(y,1));

i=1;
while t<tEnd
  Y(:,1)=y;
  for j=2:RK.stage
    MY(j-1).M=MMatrix(Reak,Species,t,Y(:,j-1));
    M=sparse(size(y,1),size(y,1));
    for k=1:j-1
      M=M+RK.a(j,k)*MY(k).M;
    end
    %Y(:,j)=(eye(size(y,1))-dt*M)\y;
    Y(:,j)=expm(dt*M)*y;
  end
  MY(RK.stage).M=MMatrix(Reak,Species,t,Y(:,RK.stage));
  M=sparse(size(y,1),size(y,1));
  for j=1:RK.stage
    M=M+RK.b(j)*MY(j).M;
  end
  %y=(eye(size(y,1))-dt*M)\y;
  y=expm(dt*M)*y;
  t=t+dt;
  i=i+1;
  TOut(i+1)=t;
  YOut(:,i+1)=y;
  dt=min(2*dt,dtMax);
  dt=min(dt,tEnd-t);
end
TOut=TOut';
YOut=YOut';
end


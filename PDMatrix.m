function [P,D]=PDMatrix(t,y)
global Bsp
switch Bsp
    case 'Brusselator'
        [P,D]=BrusselatorPD(y);
    case 'NonlinTest'
        [P,D]=NonlinTestPD(y);
end
end

function [P,D]=BrusselatorPD(y)
P=zeros(size(y,1),size(y,1));
D=zeros(size(y,1),size(y,1));

% p32(y) = d23(y) = k2y2y5, p45(y) = d54(y) = k4y5, p51(y) = d15(y) = k1y1,
% p56(y) = d65(y) = k3y2^2y6, p65(y) = d56(y) = k2y2y5,

P(3,2)=y(2)*y(5);
P(4,5)=y(5);
P(5,1)=y(1);
P(5,6)=y(5)^2*y(6);
P(6,5)=y(2)*y(5);


D(2,2)=D(2,2)+P(3,2);
D(5,5)=D(5,5)+P(4,5);
D(1,1)=D(1,1)+P(5,1);
D(6,6)=D(6,6)+P(5,6);
D(5,5)=D(5,5)+P(6,5);
end


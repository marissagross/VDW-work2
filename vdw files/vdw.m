function dy = vdw(t,y)

global eps;
global D;
global d;
%global x;
global m;
global sigma;
global gamma;
global eta;
global w;
dy = zeros(2,1);

dy(1) = y(2);

%dy(2) = ((12*eps/(m*sigma))*((sigma/y(1))^13 - (sigma/y(1))^7 - (sigma/(D-y(1)))^13 + (sigma/(D-y(1)))^7) - gamma*y(2));
dy(2) = w*(((12*y(1)*(sigma)^(12))/(y(1)^2 + d^2)^7) - ((12*y(1)*(sigma)^(6))/(y(1)^2 + d^2)^4)) +  w*(((12*y(1)*(sigma)^(12))/((y(1)^2 + (D-d)^2)^7))- ((12*y(1)*(sigma)^(6))/((y(1)^2 + (D-d)^2)^4))) - y(2)*eta;
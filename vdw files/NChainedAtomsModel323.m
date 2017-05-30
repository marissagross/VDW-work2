%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The right-hand side of the equation is defined in 
% the matlab function file 'NChainedAtoms.m'. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hl = 1;      % distance the left fixed atoms are spaced apart
hr = 1;      % distance the right fixed atoms are spaced apart
H = 0;       % offsets the right wall from the origin
Nafix = 5;   % number of atoms above and below the free atoms interact with at a time
Nafree = 3;  % number of free atoms in the chain
eta = 1;     % friction coefficient
D = 4;       % distance between walls. Here the origin is set at 0. the walls are at +- D/2
sigma = 1;   % 'happy distance'
w = 1;       
k = 1;       % spring constant
l = 1;       % natural length of the spring

t = [0 3];   % Define the time interval over which solution will be computed.  

% read this as:
x1 = 1.5;        % init x1
y1 = 1.5;      % init y1
x1p = 0;       % init x1'
y1p = 0;       % init y1'

x2 = 0; 
y2 = 0;       
x2p = 0;      
y2p = 0;      

x3 = -1.5;
y3 = -1.5;
x3p = 0;
y3p = 0;

%set initial conditions
init = [x1 y1 x2 y2 x3 y3 x1p y1p x2p y2p x3p y3p];
% Set some options used in the next command. 
options = odeset('RelTol',1e-8, 'AbsTol',1e-8);

% Numerically solve the equation.  See 'ode45' in matlab help for more info on this commmand. 
[T,Y] = ode45(@NChainedAtoms323,t,init,options,eta,D,w,sigma,hl,hr,H,k,l,Nafix,Nafree);


% Plot the solution. 

set(gca,'FontSize',24)

yy = -25:.005:25;

hold on;

for k = 1:length(yy)
    plot(-D/2, yy(k), '.b')  % plots left wall
    plot(D/2, yy(k), '.b')   % plots right wall
end
plot(Y(:,1), Y(:,2),'r') % plots first atom in red
plot(Y(:,3), Y(:,4),'k') % plots second atom in black
plot(Y(:,5), Y(:,6),'b') % plots third atom in blue

plot(Y(end,1),Y(end,2),'o'); % plot where the end up
plot(Y(end,3),Y(end,4),'o');
plot(Y(end,5),Y(end,6),'o'); 

hold off;

figure;
plot(T,Y(:,1),T,Y(:,2))
legend('x1(t)','y1(t)')

figure;
plot(T,Y(:,3),T,Y(:,4))
legend('x2(t)','y2(t)')

figure;
plot(T,Y(:,5),T,Y(:,6))
legend('x3(t)','y3(t)')
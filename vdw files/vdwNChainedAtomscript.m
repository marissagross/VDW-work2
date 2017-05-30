%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The right-hand side of the equation is defined in 
% the matlab function file 'vdw.m'. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hl = 1;      % distance the left fixed atoms are spaced apart
hr = 1;      % distance the right fixed atoms are spaced apart
H = 0;       % offsets the right wall from the origin
Nafix = 5;      % number of atoms above and below the free atoms interact with at a time
Nafree = 4;
eta = 1;     % friction coefficient
D = 4;       % distance between walls. Here the origin is set at 0. the walls are at +- D/2
sigma = 1;   % 'happy distance'
w = 1;       % not sure what the interpretation of this is

% distance from all left wall fixed atoms
% sqrt(x^2 + (ky + j*hl)^2)

% total van der waals with left wall
% summation of V(dj) from -inf < j < inf

% Use Na for number of atoms needed so we do not need an infinite number of
% them. This will depend on hl and sigma

% new summation:
% summation of V(dj) from -Na < j < Na

% ky = mod(y, hl)

% we also need the forces the two free atoms impose on each other.
% we model this using L-J 12/6 as a spring

t = [0 20];   % Define the time interval over which solution will be computed.  

% read this as:
x1 = D/2;     % init x1
y1 = 3.1;     % init y1
x1p = 0;      % init x1'
y1p = 0;      % init y1'
x2 = D/4;     % init x2
y2 = 1.9;     % init y2
x2p = 0;      % init x2'
y2p = 0;      % init y2'
x3 = 0;
y3 = 0;
x3p = 1;
y3p = 1;
x4 = 0;
y4 = 0;
x4p = 0;
y4p = 0;

%set initial conditions
init = [x1 y1 x2 y2 x3 y3 x4 y4 x1p y1p x2p y2p x3p y3p x4p y4p];
% Set some options used in the next command. 
options = odeset('RelTol',1e-8, 'AbsTol',1e-8);

% Numerically solve the equation.  See 'ode45' in matlab help for more info on this commmand. 
[T,Y] = ode45(@vdwNChainedAtoms,t,init,options,eta,D,w,sigma,hl,hr,H,Nafix, Nafree);


% Plot the solution. 

set(gca,'FontSize',24)

yy = -25:.005:25;

hold on;

for k = 1:length(yy)
    plot(-D/2, yy(k), '.b')  % plots left wall
    plot(D/2, yy(k), '.b')   % plots right wall
end
plot(Y(:,1), Y(:,2)) % plots first atom in red
plot(Y(:,5), Y(:,6)) % plots second atom in black
plot(Y(:,9), Y(:,10)) 
plot(Y(:,13), Y(:,14)) 
plot(Y(end,1),Y(end,2),'o'); % plot where the atoms end up
plot(Y(end,5),Y(end,6),'o');
plot(Y(end,9),Y(end,10),'o'); % plot where the atoms end up
plot(Y(end,13),Y(end,14),'o');

hold off;

figure;
plot(T,Y(:,1),T,Y(:,2),T,Y(:,5),T,Y(:,6),T,Y(:,9),T,Y(:,10),T,Y(:,13),T,Y(:,14))
legend('x1(t)','y1(t)','x2(t)','y2(t)','x3(t)','y3(t)','x4(t)','y4(t)')

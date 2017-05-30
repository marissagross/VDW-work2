%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The right-hand side of the equation is defined in 
% the matlab function file 'vdw.m'. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hl = 1;      % distance the left fixed atoms are spaced apart
hr = 1;      % distance the right fixed atoms are spaced apart
H = 0;       % offsets the right wall from the origin
Na = 6;      % number of atoms above and below the free atoms interact with at a time
eta = 0;     % friction coefficient
D = 1.8;       % distance between walls. Here the origin is set at 0. the walls are at +- D/2
sigma = 1;   % 'happy distance'
w = 1;       % the force of the walls
k = 10;       % spring force
l = 1;       % natural length of the spring 

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

t = [0 4];   % Define the time interval over which solution will be computed.  

hold on;
% read this as:
x1 = 0;     % init x1
y1 = 2.1;       % init y1
x1p = 0;      % init x1'
y1p = 100;      % init y1'
x2 = 0;      % init x2
y2 = 2;     % init y2
x2p = 0;      % init x2'
y2p = -100;      % init y2'

%set initial conditions
init = [x1 y1 x1p y1p x2 y2 x2p y2p];
% Set some options used in the next command. 
options = odeset('RelTol',1e-8, 'AbsTol',1e-8);

% Numerically solve the equation.  See 'ode45' in matlab help for more info on this commmand. 
[T,Y] = ode45(@vdw2chainmodel,t,init,options,eta,D,w,k,l,sigma,hl,hr,H,Na);


% Plot the solution. 

set(gca,'FontSize',24)

yy = -25:.005:25;

for k = 1:length(yy)
    plot(-D/2, yy(k), '.b')  % plots left wall
    plot(D/2, yy(k), '.b')   % plots right wall
end
plot(Y(:,1), Y(:,2),'r') % plots first atom in red
plot(Y(:,5), Y(:,6),'k') % plots second atom in black
plot(Y(end,1),Y(end,2),'o'); % plot where the end up
plot(Y(end,5),Y(end,6),'o');
figure(2)
plot(T,Y(:,1),T,Y(:,2),T,Y(:,5),T,Y(:,6))
legend('x1(t)','y1(t)','x2(t)','y2(t)')
hold off;

% provide the final distance between the two free atoms
r = sqrt( (Y(end,1)-Y(end,5))^2 + (Y(end,2)-Y(end,6))^2 )


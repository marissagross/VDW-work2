%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script numerically solves 2nd-order ODE
% 
% ((12*eps/m*sigma)*((sigma/x)^13 - (sigma/x)^7 - (sigma/D-x)^13 + (sigma/D-x)^7))
%
% and then plots the solution.
%
% The right-hand side of the equation is defined in 
% the matlab function file 'vdw.m'. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global eta; %friction coefficient
global D; %D in eqn
global w; %vdw strength
global sigma; %'Happy' distance

global hl;
global hr;
global H;
global Na;

%distance from all left wall fixed atoms
% sqrt(x^2 + (ky + j*hl)^2)

%total van der waals with left wall
%summation of V(dj) from -inf < j < inf

%Use Na for number of atoms needed so we do not need an infinite number of
%them. This will depend on hl and sigma

%new summation:
%summation of V(dj) from -Na < j < Na

%ky = mod(y, hl)

hl = 1;     % distance between atoms on left wall
hr = 1;     % distance between atoms on right wall
H = 0;      % distance of atom on right side of wall above x-axis
Na = 7;     % # of atoms up and down applying force on free atom

eta = 1;       % friction
D = 1.5;       % distance between walls
sigma = 1;     % "happy distance"
w = 1;         % vdw strength

t = [0 20];     % Define the time interval over which solution will be computed.  
                  % You may need to change the right end point to see the long-term
                  % behavior of the solution, depending on epsilon. 
h = animatedline;

hold on;
for j = 0:1:0
    deltav = .5;
    %d = d+.5*j;
    init = [0.25 .75 0 100];   % init x position, init y position, init x velocity, init y velocity


% Set some options used in the next command.  Do not worry about this for now.
options = odeset('RelTol',1e-8, 'AbsTol',1e-8);

% Numerically solve the equation.  See 'ode45' in matlab help for more info on this commmand. 
[T,Y] = ode45(@vdwXY_infAtoms,t,init,options);


% Plot the solution. 

% position
%subplot(2,1,1)
%xlim([0,D])
%ylim([-3,3])
set(gca,'FontSize',24)

yy = -3:.005:3;

for k = 1:length(yy)
    plot(-D/2, yy(k), '.b')  % plots left wall
    plot(D/2, yy(k), '.b')   % plots right wall
end

 plot(Y(:,1), Y(:,2)) % plot solution

%for i = 1:length(Y(:,1))/2  %loop plots the animated line of the solution
 %   addpoints(h, Y(2*i,1), Y(2*i,2));
 %   addpoints(h, Y((2*i)+1,1), Y((2*i)+1,2));
 %   drawnow
%end

plot(Y(end, 1), Y(end, 2), 'o')

% legend('Position','Velocity')
end
hold off;
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
global h1; %y distance of atom 1 from x-axis
global h2; %y distance of atom 2 from x-axis
global h3; %y distance of atom 3 from x-axis
global h4; %y distance of atom 4 from x-axis
global h5; %y distance of atom 5 from x-axis
global h6; %y distance of atom 6 from x-axis

eta = 1;     % friction
D = 2;        %distance between walls
sigma = 1;
w = 1;

h1 = -1;
h2 = 0;
h3 = 1;
h4 = -1;
h5 = 0;
h6 = 1;

t = [0 50];     % Define the time interval over which solution will be computed.  
                  % You may need to change the right end point to see the long-term
                  % behavior of the solution, depending on epsilon. 
h = animatedline;

hold on;
for j = 0:1:0
    deltav = .5;
    %d = d+.5*j;
    init = [.25 .75 -1 0];   % init x position, init y position, init x velocity, init x velocity


% Set some options used in the next command.  Do not worry about this for now.
options = odeset('RelTol',1e-8, 'AbsTol',1e-8);

% Numerically solve the equation.  See 'ode45' in matlab help for more info on this commmand. 
[T,Y] = ode45(@vdwXY_6atoms,t,init,options);


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
plot(-D/2, h1, '.r', 'markersize',20) % plot fixed atom 1
plot(-D/2, h2, '.r', 'markersize',20) % plot fixed atom 2
plot(-D/2, h3, '.r', 'markersize',20) % plot fixed atom 3 
plot(D/2, h4, '.r', 'markersize',20) % plot fixed atom 4
plot(D/2, h5, '.r', 'markersize',20) % plot fixed atom 5
plot(D/2, h6, '.r', 'markersize',20) % plot fixed atom 6 

% plot(Y(:,1), Y(:,2)) % plot solution

for i = 1:length(Y(:,1))/2  %loop plots the animated line of the solution
    addpoints(h, Y(2*i,1), Y(2*i,2));
    addpoints(h, Y((2*i)+1,1), Y((2*i)+1,2));
    drawnow
end

plot(Y(end, 1), Y(end, 2), 'o')

% legend('Position','Velocity')
end
hold off;
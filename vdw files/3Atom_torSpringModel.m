%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This solves the problem for a chain of atoms placed between two walls.
% This code allows you to choose how many atoms are in the chain freely.
% This code assumes mass = 0, so it reduces to a first order problem.
% For each atom, you will need 2 IC.
% The right-hand side of the equation is defined in 
% the matlab function file 'NAtom_torSpring.m'. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paramaters
tic;
hl = 1;      % distance the left fixed atoms are spaced apart
hr = 1;      % distance the right fixed atoms are spaced apart
H = 0;       % offsets the right wall from the origin
Nafix = 10;  % number of atoms above and below the free atoms interact with at a time
Nafree = 7;  % number of free atoms in the chain
eta = 10;    % friction coefficient
D = 4;       % distance between walls. Here the origin is set at 0. the walls are at +- D/2
sigma = 1;   % 'happy distance' - vdw potential is zero at sigma
w = 1;       % strength of vdw interaction
k = 1;       % spring constant for chain of atoms
l = 1;       % natural length of the spring for chain of atoms
mu = 1;      % torsional spring coefficient

t = [0 120];   % Define the time interval over which solution will be computed.  

% set initial conditions
    init = zeros(1,2*Nafree);
    % position loop
    for i = 1:Nafree
        init((2*i) - 1) = 0.5*(-1)^i;    % x initial position
        init(2*i) = 1.5*i;    % y initial positions
    end

% Set some options used in the next command. 
options = odeset('RelTol',1e-8, 'AbsTol',1e-8);

% Numerically solve the equation.  See 'ode45' in matlab help for more info on this commmand. 
[T,Y] = ode45(@NAtom_torSpring,t,init,options,eta,mu,D,w,sigma,hl,hr,H,k,l,Nafix,Nafree);

% Plot the solution. 

set(gca,'FontSize',24)

yyl = -25:hl:25;
yyr = -25:hr:25;

hold on;
for k = 1:length(yyl)
    plot(-D/2, yyl(k), '.b')  % plots left wall
    plot(D/2, yyr(k)+H, '.b')   % plots right wall
end

% parametric plots
for i = 1:Nafree
    plot(Y(:,(2*i)-1),Y(:,2*i)) % plots its path
    plot(Y(end,(2*i)-1),Y(end,2*i),'o') % plots where they end up
end
hold off;

% spit out the length of the springs to compare with natural length
r = zeros(1,Nafree-1);
for i = 1:Nafree-1
    r(i) = norm([Y(end,2*i-1),Y(end,2*i)]-[Y(end,2*i+1),Y(end,2*i+2)]);
end
axis equal;
disp(r);
toc; 
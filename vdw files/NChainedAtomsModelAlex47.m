%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This solves the problem for a chain of atoms placed between two walls.
% This code allows you to choose how many atoms are in the chain freely.
% For each atom, you will need 4 IC.
% The right-hand side of the equation is defined in 
% the matlab function file 'NChainedAtoms.m'. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paramaters
tic;
hl = 2;      % distance the left fixed atoms are spaced apart
hr = 2;      % distance the right fixed atoms are spaced apart
H = 0;       % offsets the right wall from the origin
Nafix = 15;   % number of atoms above and below the free atoms interact with at a time
Nafree = 4;  % number of free atoms in the chain
eta = 5;    % friction coefficient
D = 4;       % distance between walls. Here the origin is set at 0. the walls are at +- D/2
sigma = 1;   % 'happy distance' - vdw potential is zero at sigma
w = 1;       % strength of vdw interaction
k = 1;      % spring constant for chain of atoms
l = 1.5;       % natural length of the spring for chain of atoms

t = [0 15];   % Define the time interval over which solution will be computed.  
 
% set initial conditions
init = zeros(1,4*Nafree);

% initial position loop
for i = 1:Nafree
    init((2*i)-1) = 0; % x initial position
    init(2*i) = i;  % y initial position
end

%initial velocity loop
for i = 1:Nafree
    init(2*Nafree+(2*i)-1) = 1; % x initial position
    init(2*Nafree+(2*i)) = 0; % y initial position
end

% Set some options used in the next command. 
options = odeset('RelTol',1e-8, 'AbsTol',1e-8);

% Numerically solve the equation.  See 'ode45' in matlab help for more info on this commmand. 
[T,Y] = ode45(@NChainedAtomsAkex47,t,init,options,eta,D,w,sigma,hl,hr,H,k,l,Nafix,Nafree);

% Plot the solution. 

set(gca,'FontSize',24)

yy = -25:.005:25;

hold on;

for k = 1:length(yy)
    plot(-D/2, yy(k), '.b')  % plots left wall
    plot(D/2, yy(k), '.b')   % plots right wall
end

% parametric plots
for i = 1:Nafree
    plot(Y(:,2*i-1),Y(:,2*i)) % plots its path
    plot(Y(end,2*i-1),Y(end,2*i),'o') % plots where they end up
end
hold off;

% spit out their locations from each other in order to gauge if they are at
% equilibrium position for the final time step
% this tells you the distance each is from their neighbors
r = zeros(1,Nafree-1);
for i = 1:Nafree-1
    r(i) = norm([Y(end,2*i-1),Y(end,2*i)]-[Y(end,2*i+1),Y(end,2*i+2)]);
end
disp(r);
toc; 
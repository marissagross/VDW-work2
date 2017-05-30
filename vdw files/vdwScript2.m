% x'' = 12*eps/m*sigma)*((sigma/x)^13 - (sigma/x)^7 - (sigma/D-x)^13 + (sigma/D-x)^7))

% (1) illustrates the use of the ode solvers

eps = 1;          % Set the value of epsilon.
D = 1;
m = 1;
sigma = D/2;
x = sigma;

y0 = [0 D];
tol = 1e-6;
tspan = [0 100];

options = odeset('AbsTol',tol,'RelTol',tol);

[t,y] = ode45(@vdw,tspan,y0,options);


figure(1)
plot(t,y) %,legend('ode45')


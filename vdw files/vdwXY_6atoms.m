function dy = vdw(t,y)

global D;
global w;
global sigma;
global eta;
global h1;
global h2;
global h3;
global h4;
global h5;
global h6;

dy = zeros(4,1);

%%

r1 = [y(1), y(2)] - [-D/2, h1]; % vector from fixed atom 1 - left wall
r2 = [y(1), y(2)] - [-D/2, h2]; % vector from fixed atom 2 - left wall
r3 = [y(1), y(2)] - [-D/2, h3]; % vector from fixed atom 3 - left wall
r4 = [y(1), y(2)] - [D/2, h4]; % vector from fixed atom 4 - right wall
r5 = [y(1), y(2)] - [D/2, h5]; % vector from fixed atom 5 - right wall
r6 = [y(1), y(2)] - [D/2, h6]; % vector from fixed atom 6 - right wall

r1n = norm(r1);
r2n = norm(r2);
r3n = norm(r3);
r4n = norm(r4);
r5n = norm(r5);
r6n = norm(r6);

% force exerted on free atom by fixed atom 1
f1 = 12*w*(sigma^12/r1n.^13 - sigma^6/r1n.^7)*r1/r1n;

% force exerted on free atom by fixed atom 2
f2 = 12*w*(sigma^12/r2n.^13 - sigma^6/r2n.^7)*r2/r2n;

% force exerted on free atom by fixed atom 3
f3 = 12*w*(sigma^12/r3n.^13 - sigma^6/r3n.^7)*r3/r3n;

% force exerted on free atom by fixed atom 4
f4 = 12*w*(sigma^12/r4n.^13 - sigma^6/r4n.^7)*r4/r4n;

% force exerted on free atom by fixed atom 5
f5 = 12*w*(sigma^12/r5n.^13 - sigma^6/r5n.^7)*r5/r5n;

% force exerted on free atom by fixed atom 6
f6 = 12*w*(sigma^12/r6n.^13 - sigma^6/r6n.^7)*r6/r6n;


dy(1) = y(3);

dy(2) = y(4);

dy(3) = f1(1) + f2(1) + f3(1) + f4(1) + f5(1) + f6(1) - eta*y(3);

dy(4) = f1(2) + f2(2) + f3(2) + f4(2) + f5(2) + f6(2) - eta*y(4);

end

Bending Forces

sin is approximately x when x is close to zero
tan also is approximately x when x is close to zero

torisonal spring

deltaTheta = theta(i-1) - theta(i)
deltaTheta wants to be zero

bending energy = (mu/2)*(deltaTheta)^2

= 2*mu*tan^2(deltaTheta/2)
= 2*mu*(sin^2(deltaTheat/2)/cos^2(deltaTheta))

= 2*mu*((1-cos(deltaTheta))/(1+cos(deltaTheta)))

= 2*mu*((norm(r(i)-r(i-1))*norm(r(i+1)-r(i))-(r(i)-r(i-1))o(r(i+1)-r(i)))/(norm(r(i)-r(i-1))*norm(r(i+1)-r(i))+(r(i)-r(i-1))o(r(i+1)-r(i))))

1) compute derivatives of bending energy to get bending forces
2) write down the correct equations for three atoms with one torsional spring
3) run simulations of 2) and show results

a = r(i) - r(i-1);
b = r(i+1) - r(i);

aN = norm(a);
bN = norm(b);

top = (aN*bN)-dot(a,b);
bottom = (aN*bN)+dot(a,b);

topP = (a(1))*(bN/aN) + (-b(1))*(aN/bN) - (b(1) - a(1));
bottomP = topP + 2*(b(1) - a(1));

torSpr(i,1) = 2*mu*(topP*bottom - top*bottomP)/(bottom^2);

topP = (a(2))*(bN/aN) + (-b(2))*(aN/bN) - (b(2) - a(2));
bottomP = topP + 2*(b(2) - a(2));

torSpr(i,2) = 2*mu*(topP*bottom - top*bottomP)/(bottom^2);

topP = (-a(1))*(bN/aN) + b(1);
bottomP = topP - 2*b(1);

torSpr(i-1,1) = 2*mu*(topP*bottom - top*bottomP)/(bottom^2);

topP = (-a(2))*(bN/aN) + b(2);
bottomP = topP - 2*b(2);

torSpr(i-1,2) = 2*mu*(topP*bottom - top*bottomP)/(bottom^2);

topP = (b(1))*(bN/aN) - a(1);
bottomP = topP + 2*a(1);

torSpr(i+1,1) = 2*mu*(topP*bottom - top*bottomP)/(bottom^2);

topP = (b(2))*(bN/aN) - a(2);
bottomP = topP + 2*a(2);

torSpr(i+1,2) = 2*mu*(topP*bottom - top*bottomP)/(bottom^2);
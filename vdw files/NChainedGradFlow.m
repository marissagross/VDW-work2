function dy = NChainedGradFlow(t,y,eta,D,w,sigma,hl,hr,H,k,l,Nafix,Nafree)
% the ordering this will spit out is:
% Y = x1,y1,x2,y2,...,xn,yn
% we will have 2 equations per free atom
dy = zeros(2*Nafree,1);

% preallocate memory for speed
Fspring = zeros(Nafree, 2);
F = zeros(Nafree, 2);
kl = zeros(1, Nafree);
kr = zeros(1, Nafree);

%%
% initialize vectors pointing from walls to each free atom
for i = 1:(Nafree)
    kl(i) = mod(y(2*i), hl);
    kr(i) = mod(y(2*i) - H, hr);
end

% forces from the walls
for i = 1:Nafree
    for j = 0:2*Nafix
        %creates the forces on each free atom by both fixed walls
        rl = [y(2*i-1), 0] - [-D/2, (j - Nafix)*hl - kl(i)];
        rln = norm(rl);
        fl = 12*w*(sigma^12/rln.^13 - sigma^6/rln.^7)*rl/rln;

        rr = [y(2*i-1), 0] - [D/2, (j - Nafix)*hl - kr(i)];
        rrn = norm(rr);
        fr = 12*w*(sigma^12/rrn.^13 - sigma^6/rrn.^7)*rr/rrn;
        
        %store these in arrays to be saved for newton's second law
        F(i,:) = F(i,:) + fl + fr; % so it stores all of the left forces for free atom i in one spot 
    end
end

% create forces that the free atoms exert on each other.
% the free atoms only interact with the one above and below it
% We model these by Hooke Springs
% Fpring is indexed by x,y then free atom i then 1 for force coming from
% below, 2 for force coming from above
for i = 1:(Nafree - 1)
    Rspring = [y(2*i-1),y(2*i)]-[y(2*i+1),y(2*i+2)]; % 'spring' vector from atom above
    Fspring(i,:) = k*(l-norm(Rspring))*(Rspring)/norm(Rspring);
end

for i = Nafree:-1:2
    Fspring(i,:) = Fspring(i,:) - Fspring((i-1),:);
end

% create the diff eq's:
for i = 1:Nafree
        dy(2*i-1) = dy(2*i-1) + F(i,1) + Fspring(i,1);
        dy(2*i) = dy(2*i) + F(i,2) + Fspring(i,2);
end

% finally throw in friction
for i = 1:Nafree
    dy(2*i-1) = dy(2*i-1)/eta;
    dy(2*i) = dy(2*i)/eta;
end

end
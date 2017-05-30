function dy = NChainedAtoms323(t,y,eta,D,w,sigma,hl,hr,H,k,l,Nafix,Nafree)
% the ordering this will spit out is:
% Y = x1,y1,x2,y2,...,xn,yn,x1',y1',...xn',yn'
% we will have 4 equations per free atom
dy = zeros(4*Nafree,1);

% preallocate memory for speed
Fspring = zeros(2, Nafree, 2);
FL = zeros(2, Nafree, 2*Nafix + 1);
FR = zeros(2, Nafree, 2*Nafix + 1);
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
        FL(:,i,j+1) = fl; % so it stores all of the left forces for free atom i in one spot 
        FR(:,i,j+1) = fr; % these are multidimensional arrays
    end
end

% create forces that the free atoms exert on each other.
% the free atoms only interact with the one above and below it
% We model these by Hookeen Springs
% Fpring is indexed by x,y then free atom i then 1 for force coming from
% below, 2 for force coming from above
for i = 1:Nafree
        if i==1
            RspringA = [y(1),y(2)]-[y(3),y(4)]; % 'spring' vector from atom above
            Fspring(:,i,1) = [0,0];
            Fspring(:,i,2) = k*(l-norm(RspringA))*(RspringA)/norm(RspringA);
        elseif i==Nafree
            RspringB = [y(2*Nafree-1),y(2*Nafree)]-[y(2*Nafree-3),y(2*Nafree-2)]; % 'spring' vector from atom below
            Fspring(:,i,1) = k*(norm(RspringB)-l)*(RspringB)/norm(RspringB);
            Fspring(:,i,2) = [0,0];
        else
            RspringA = [y(2*i-1),y(2*i)]-[y(2*i+1),y(2*i+2)]; % 'spring' vector from atom above
            RspringB = [y(2*i-1),y(2*i)]-[y(2*i-3),y(2*i-2)]; % 'spring' vector from atom below
            Fspring(:,i,1) = k*(norm(RspringB)-l)*(RspringB)/norm(RspringB);
            Fspring(:,i,2) = k*(l-norm(RspringA))*(RspringA)/norm(RspringA);
        end
end
% create the first diff eqs for all the positions
for i = 1:Nafree
    dy(2*i-1) = y(2*Nafree+(2*i-1)); % the check is dy(1) = y(2Nafree+1) which is the x1,1' = x1,2
    dy(2*i) = y(2*Nafree+(2*i));
end

% now we work on making the next set for the velocities
% sum up the forces from the fixed atoms and throw them in to the eqn
for i = 1:Nafree
    for p = 0:2*Nafix
    
        dy(2*Nafree+(2*i-1)) = dy(2*Nafree+(2*i-1)) + FL(1,1,p+1) + FR(1,1,p+1);
        dy(2*Nafree+(2*i)) = dy(2*Nafree+(2*i)) + FL(1,2,p+1) + FR(1,2,p+1);
    end
end

% now we must also use their forces on each other:
for i = 1:Nafree
        dy(2*Nafree+(2*i-1)) = dy(2*Nafree+(2*i-1)) + Fspring(1,i,1) + Fspring(1,i,2);
        dy(2*Nafree+(2*i)) = dy(2*Nafree+(2*i)) + Fspring(2,i,1) + Fspring(2,i,2);
end

% finally throw in friction
for i = 1:Nafree
    dy(2*Nafree+(2*i-1)) = dy(2*Nafree+(2*i-1)) -eta*y(2*Nafree+(2*i-1));
    dy(2*Nafree+(2*i)) = dy(2*Nafree+(2*i)) - eta*y(2*Nafree+(2*i));
end

end
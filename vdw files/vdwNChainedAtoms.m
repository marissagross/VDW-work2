function dy = vdwNChainedAtoms(t,y,eta,D,w,sigma,hl,hr,H,Nafix,Nafree)
% the ordering this will spit out is:
% Y = x1,y1,x2,y2,...,xn,yn,x1',y1',...xn',yn'
% we will have 4 equations per free atom
dy = zeros(4*Nafree,1);
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
for k = 1:Nafree
    for j = 0:2*Nafix
        %creates the forces on free atom 1 by both fixed walls
        rl = [y(2*k-1), 0] - [-D/2, (j - Nafix)*hl - kl(k)];
        rln = norm(rl);
        fl = 12*w*(sigma^12/rln.^13 - sigma^6/rln.^7)*rl/rln;

        rr = [y(2*k-1), 0] - [D/2, (j - Nafix)*hl - kr(k)];
        rrn = norm(rr);
        fr = 12*w*(sigma^12/rrn.^13 - sigma^6/rrn.^7)*rr/rrn;
        
        %store these in vectors to be saved for newton's second law
        FL(:,k,j+1) = fl; % so it stores all of the left forces for free atom k in one spot 
        FR(:,k,j+1) = fr; % these are multidimensional arrays
    end
end

% create forces that the free atoms exert on each other.
% the free atoms only interact with the one above and below it
for i = 1:Nafree
        if i==1
            RspringA = [y(2*i-1),y(2*i)]-[y(2*i+1),y(2*i+2)]; % 'spring' vector from atom above
            Fspring(:,i,1) = [0,0];
            Fspring(:,i,2) = 12*w*(sigma^12/norm(RspringA).^13 - sigma^6/norm(RspringA).^7)*RspringA/norm(RspringA);
        elseif i==Nafree
            RspringB = [y(2*i-1),y(2*i)]-[y(2*i-3),y(2*i-2)]; % 'spring' vector from atom below
            Fspring(:,i,1) = 12*w*(sigma^12/norm(RspringB).^13 - sigma^6/norm(RspringB).^7)*RspringB/norm(RspringB);
            Fspring(:,i,2) = [0,0];
        else
            RspringA = [y(2*i-1),y(2*i)]-[y(2*i+1),y(2*i+2)]; % 'spring' vector from atom above
            RspringB = [y(2*i-1),y(2*i)]-[y(2*i-3),y(2*i-2)]; % 'spring' vector from atom below
            Fspring(:,i,1) = 12*w*(sigma^12/norm(RspringB).^13 - sigma^6/norm(RspringB).^7)*RspringB/norm(RspringB);
            Fspring(:,i,2) = 12*w*(sigma^12/norm(RspringA).^13 - sigma^6/norm(RspringA).^7)*RspringA/norm(RspringA);
        end
end
% create the first diff eqs for all the positions
for i = 1:Nafree
    dy(2*i-1) = y(2*i-1);
    dy(2*i) = y(2*i);
end

% now we work on making the next set for the velocities
% sum up the forces from the fixed atoms and throw them in to the eqn
for i = (Nafree+1):2*Nafree
    for p = 0:2*Nafix
    
        dy(2*i-1) = dy(2*i-1) + FL(1,1,p+1) + FR(1,1,p+1);
        dy(2*i) = dy(2*i) + FL(1,2,p+1) + FR(1,2,p+1);
    end
end

% now we must also use their forces on each other:
for i = 1:Nafree
        dy(2*i-1+Nafree) = dy(2*i-1+Nafree) + Fspring(1,i,1) + Fspring(1,i,2) - eta*y(2*i-1+Nafree); % added friction too
        dy(2*i+Nafree) = dy(2*i+Nafree) + Fspring(2,i,1) + Fspring(2,i,2) - eta*y(2*i+Nafree);    
end
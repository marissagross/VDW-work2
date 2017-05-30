function dy = vdwlongchainmodel(t,y,eta,D,w,sigma,hl,hr,H,Nafix,Nafree)
% the ordering this will spit out is:
% Y = x1,y1,x2,y2,...,xn,yn,x1',y1',...xn',yn'
% we will have 4 equations per free atom
dy = zeros(4*Nafree,1);
FL = zeros(2, Nafree, 2*Nafix + 1);
FR = zeros(2, Nafree, 2*Nafix + 1);

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
        rl(:, j + 1) = [y(2*k-1), 0] - [-D/2, (j - Nafix)*hl - kl(k)];
        rln(j + 1) = norm(rl(:,j + 1));
        fl(:,j + 1) = 12*w*(sigma^12/rln(j+1).^13 - sigma^6/rln(j+1).^7)*rl(:,j+1)/rln(j+1);

        rr(:, j + 1) = [y(2*k-1), 0] - [D/2, (j - Nafix)*hl - kr(k)];
        rrn(j + 1) = norm(rr(:,j + 1));
        fr(:,j + 1) = 12*w*(sigma^12/rrn(j+1).^13 - sigma^6/rrn(j+1).^7)*rr(:,j+1)/rrn(j+1);
        %store these in vectors to be saved for newton's second law
        FL(:,k,j) = fl; % so it stores all of the left forces for free atom k in one spot 
        FR(:,k,j) = fr; % these are multidimensional arrays
    end
end

% create forces that the free atoms exert on each other.
% the free atoms only interact with the one above and below it
for i = 1:Nafree
% create the first diff eqs for all the positions
for i = 1:2*Nafree
    dy(2*i-1) = y((2*Na)+(2*i-1));
    dy(2*i) = y((2*Na)+(2*i));
end
end

% now we work on making the next set for the velocities
% sum up the forces from the fixed atoms and throw them in to the eqn
for i = (2*Nafree+1):4*Nafree
    for p = 0:2*Na
    
        dy((2*Nafree)+(2*i-1)) = dy((2*Nafree)+(2*i-1)) + FL(1,1,p+1) + FR(1,1,p+1);
        dy((2*Nafree)+(2*i)) = dy((2*Nafree)+(2*i)) + FL(1,2,p+1) + FR(1,2,p+1);
    end
end

% now we must also use their forces on each other:
for i = (2*Nafree+1):4*Nafree
    
end

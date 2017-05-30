function dy = vdw2chainmodel(t,y,eta,D,w,sigma,hl,hr,H,Na)

dy = zeros(12,1); %if we have three connected atoms we will have 12 diff eq's

%%
u1 = [y(1),y(2)];       %represents x1, y1 i.e. position for free atom 1
u2 = [y(5),y(6)];       %represents x2, y2, position for free atom 2
u3 = [y(9),y(10)];      %represents x3, y3, position for free atom 3
kl1 = mod(y(2), hl);    %for free atom 1
kr1 = mod(y(2) - H, hr);
kl2 = mod(y(6), hl);    %for free atom 2
kr2 = mod(y(6) - H, hr);
kl3 = mod(y(10), hl);    %for free atom 3
kr3 = mod(y(10) - H, hr);

%create f on 1 by 2
F12 = 12*w*(sigma^12/norm((u1-u2))^13 - sigma^6/(norm(u1-u2))^7)*(u1-u2)/norm(u1-u2);
%create f on 1 by 3
F13 = 12*w*(sigma^12/norm((u1-u3))^13 - sigma^6/(norm(u1-u3))^7)*(u1-u3)/norm(u1-u3);
%create f on 2 by 1
F21 = 12*w*(sigma^12/norm((u2-u1))^13 - sigma^6/(norm(u2-u1))^7)*(u2-u1)/norm(u2-u1);
%create f on 2 by 3
F23 = 12*w*(sigma^12/norm((u2-u3))^13 - sigma^6/(norm(u2-u3))^7)*(u2-u3)/norm(u2-u3);
%create f on 3 by 1
F31 = 12*w*(sigma^12/norm((u3-u1))^13 - sigma^6/(norm(u3-u1))^7)*(u3-u1)/norm(u3-u1);
%create f on 3 by 2
F32 = 12*w*(sigma^12/norm((u3-u2))^13 - sigma^6/(norm(u3-u2))^7)*(u3-u2)/norm(u3-u2);


%forces from the walls
for j = 0:2*Na
    %creates the forces on free atom 1 by both fixed walls
    rl1(:, j + 1) = [y(1), 0] - [-D/2, (j - Na)*hl - kl1];
    rln1(j + 1) = norm(rl1(:,j + 1));
    fl1(:,j + 1) = 12*w*(sigma^12/rln1(j+1).^13 - sigma^6/rln1(j+1).^7)*rl1(:,j+1)/rln1(j+1);

    rr1(:, j + 1) = [y(1), 0] - [D/2, (j - Na)*hl - kr1];
    rrn1(j + 1) = norm(rr1(:,j + 1));
    fr1(:,j + 1) = 12*w*(sigma^12/rrn1(j+1).^13 - sigma^6/rrn1(j+1).^7)*rr1(:,j+1)/rrn1(j+1);
    
    %create forces on free atom 2 by both fixed walls
    rl2(:, j + 1) = [y(5), 0] - [-D/2, (j - Na)*hl - kl2];
    rln2(j + 1) = norm(rl2(:,j + 1));
    fl2(:,j + 1) = 12*w*(sigma^12/rln2(j+1).^13 - sigma^6/rln2(j+1).^7)*rl2(:,j+1)/rln2(j+1);

    rr2(:, j + 1) = [y(5), 0] - [D/2, (j - Na)*hl - kr2];
    rrn2(j + 1) = norm(rr2(:,j + 1));
    fr2(:,j + 1) = 12*w*(sigma^12/rrn2(j+1).^13 - sigma^6/rrn2(j+1).^7)*rr2(:,j+1)/rrn2(j+1);

    %create forces on free atom 3 by both fixed walls
    rl3(:, j + 1) = [y(9), 0] - [-D/2, (j - Na)*hl - kl3];
    rln3(j + 1) = norm(rl3(:,j + 1));
    fl3(:,j + 1) = 12*w*(sigma^12/rln3(j+1).^13 - sigma^6/rln3(j+1).^7)*rl3(:,j+1)/rln3(j+1);

    rr3(:, j + 1) = [y(9), 0] - [D/2, (j - Na)*hl - kr3];
    rrn3(j + 1) = norm(rr3(:,j + 1));
    fr3(:,j + 1) = 12*w*(sigma^12/rrn3(j+1).^13 - sigma^6/rrn3(j+1).^7)*rr3(:,j+1)/rrn3(j+1);

end


dy(1) = y(3);   

dy(2) = y(4);

dy(5) = y(7);

dy(6) = y(8);

dy(9) = y(11);

dy(10) = y(12);

% sum up the forces from the fixed atoms and throw them in to the eqn
for p = 0:2*Na
    
    dy(3) = dy(3) + fl1(1,p+1) + fr1(1,p+1);
    dy(4) = dy(4) + fl1(2,p+1) + fr1(2,p+1);
    
    dy(7) = dy(7) + fl2(1,p+1) + fr2(1,p+1);
    dy(8) = dy(8) + fl2(2,p+1) + fr2(2,p+1);
    
    dy(11) = dy(11) + fl3(1,p+1) + fr3(1,p+1);
    dy(12) = dy(12) + fl3(2,p+1) + fr3(2,p+1);

end

dy(3) = dy(3) + F12(1) + F13(1) - eta*y(3);
dy(4) = dy(4) + F12(2) + F13(2) - eta*y(4);

dy(7) = dy(7) + F21(1) + F23(1) - eta*y(7);
dy(8) = dy(8) + F21(2) + F23(2) - eta*y(8);

dy(11) = dy(11) + F31(1) + F32(1) - eta*y(11);
dy(12) = dy(12) + F31(2) + F32(2) - eta*y(12);



end

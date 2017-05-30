function dy = vdw(t,y)

global D;
global w;
global sigma;
global eta;
global hl;
global hr;
global Na;
global H;


dy = zeros(4,1);
rl = zeros(2, 2*Na + 1);
rln = zeros(1, 2*Na + 1);
fl = zeros(2, 2*Na + 1);

%%

kl = mod(y(2), hl);
kr = mod(y(2) - H, hr);

for j = 0:2*Na
    
    rl(:, j + 1) = [y(1), 0] - [-D/2, (j - Na)*hl + kl];
    rln(j + 1) = norm(rl(:,j + 1));
    fl(:,j + 1) = 12*w*(sigma^12/rln(j+1).^13 - sigma^6/rln(j+1).^7)*rl(:,j+1)/rln(j+1);

    rr(:, j + 1) = [y(1), 0] - [D/2, (j - Na)*hl + kr];
    rrn(j + 1) = norm(rr(:,j + 1));
    fr(:,j + 1) = 12*w*(sigma^12/rrn(j+1).^13 - sigma^6/rrn(j+1).^7)*rr(:,j+1)/rrn(j+1);
    
end


dy(1) = y(3);

dy(2) = y(4);

for p = 0:2*Na
    
    dy(3) = dy(3) + fl(1,p+1) + fr(1,p+1);
    dy(4) = dy(4) + fl(2,p+1) + fr(2,p+1);

end

dy(3) = dy(3) - eta*y(3);
dy(4) = dy(4) - eta*y(4);

end

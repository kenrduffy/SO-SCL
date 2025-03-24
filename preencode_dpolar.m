function u = preencode_dpolar(msg, frz, dCons)
n = length(frz);
assert(length(msg)+sum(frz) == n)

u = zeros(n,1);
u(~frz) = msg;

for i = 1:size(dCons,1)
    u(dCons(i,1)) = mod(u(dCons(i,1))+u(dCons(i,2)), 2);
end
end
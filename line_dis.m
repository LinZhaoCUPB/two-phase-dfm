function y = line_dis(p, vertices, v1, v2, v3)
% Distance function of a complex fracture networks
s = inf;
for i = 1:size(vertices,1)/2
    intx = 2*i-1:2*i;
    f    = vertices(intx,:);
    s    = min(s, v1 + v2 * abs(dpoly(p, f)) );
end
s = min(s, v3);
y = s;
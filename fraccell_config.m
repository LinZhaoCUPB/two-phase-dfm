function [fcp, constraints, tag] = fraccell_config(box, vertices, space)
% Generates fracture cell points, corresponding constraints and maps
xx = linspace(box(1,1), box(2,1), 10);
yy = linspace(box(1,2), box(2,2), 10);
[X, Y] = meshgrid(xx, yy);
p = [X(:), Y(:)];

constraints = reshape(1 : size(vertices,1),2,[])' ;
constraints = [constraints, (1 : size(constraints,1))'];

args = struct('precision',1e-5);
[vertices, constraints] = removeFractureIntersections(vertices, constraints, box, args);

numOrdPt = size(p,1);
p  = [p ; vertices];
constraints(:,1:2) = constraints(:,1:2) + numOrdPt;

[p, constraints] = partition_edges(p, constraints, space, box, args);

fcp = p(numOrdPt+1:end,:);

tag = constraints(:,3);

constraints = constraints(:,1:2);
constraints = constraints - numOrdPt;
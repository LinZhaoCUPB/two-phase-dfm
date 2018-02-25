function operators = operators_dfm(G, rock)
% Sets up operators for the DFM
%
% PARAMETERS:
% G         - grid structure of dfm
% rock      - rock properties
%
%
% RETURNS:
%   operators - Operators structure. See `setupOperatorsTPFA`.
%
% Written by Lin Zhao, CUPB, CHINA
%
% SEE ALSO:
%   `setupOperatorsTPFA`, `operators_radialGrid`

T_m = computeTrans_DFM(G, rock, 'hybrid', true);

[G, T_f] = computeHybridTrans(G, T_m);

N = G.faces.neighbors;
intInx = all(N ~= 0, 2);

cf = G.cells.faces(:,1);
nf = G.faces.num;

T_m  = 1./ accumarray(cf, 1./T_m, [nf, 1]);
T_m  = T_m(intInx);


N_m  = N(intInx, :);


normcon = ~isnan(T_m);
N_m = N_m(normcon,:);
T_m = T_m(normcon);


N_f  = G.cells.neighbors;

T    = [T_m; T_f];
operators.T  = T;

pv = poreVolume(G, rock);
operators.pv = pv;

% C - (transpose) divergence matrix
N  = [N_m; N_f];            
nf = size(N,1);
nc = G.cells.num;
C  = sparse( [(1:nf)'; (1:nf)'], N, ones(nf,1)*[1 -1], nf, nc);
operators.C = C;

operators.Grad = @(x)-C * x;
operators.Div  = @(x) C'* x;

% faceAvg - as multiplication with matrix
M  = sparse((1:nf)'*[1 1], N, .5*ones(nf,2), nf, nc);
operators.M = M;
operators.faceAvg = @(x) M * x;

% faceUpstr - as multiplication with matrix
upw = @(flag, x)faceUpstr(flag, x, N, [nf, nc]);
operators.faceUpstr = upw;

operators.splitFaceCellValue = @(operators, flag, x) splitFaceCellValue(operators, flag, x, [nf, nc]);

% Include neighbor relations
operators.N = N;
operators.internalConn = intInx;


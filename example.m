%% Two-phase DFM (Discrete Fracture Model) with High-Quality Triangular Grid
clear
mrstModule add ad-core ad-blackoil ad-props mrst-gui dfm
close all
addpath('distmesh')
%% Define reservoir boundary
len  = 100; hei = 100;
% Reservoir boundary box 
box  = [0,0; len, hei];
%% Fracture configurations
% Fracture vertices: [x1,    y1;   x2,    y2; .....];
%                     endpoint1,   endpoint2
vertices = ...
    [15 ,85; 85, 70;...
     15,  5; 80, 20;...
     75,  5; 25, 95;...
     20, 25; 80, 65;...
     70, 85; 50,  5;...
     20, 40; 30, 70];
 
% Length of fracture cell (triangle edge)
space = 2;
% Generate fracture cell points, corresponding constraints and maps
[fcp, constraints, tags] = fraccell_config(box, vertices, space);
%% Mesh generation by DISTMESH
% Initial edge length
h0 = 2;

%  Define distance function of the rectangular boundary
fd = @(p)drectangle(p,0,len,0,hei);

% Define scaled edge length function fh(p):
% edge length = min(v3, v1 + v2 * distance to nearest fracture)
v1 = 2;
v2 = 0.1;
v3 = 40;
fh = @(p)line_dis(p, vertices, v1, v2, v3);

% Maximum number of iterations
maxIter = 300;

% Fix points, including four reservoir boundary points and fracure cell
% points
bd_p = [0,0; 0,hei; len,0; len,hei];
fixp = [fcp; bd_p];

% Generate points of a quality triangle mesh
p  = distmesh_2d(fd, fh, h0, box, maxIter, fixp);
%% Define triangle grid by delaunayTriangulation
delTri = delaunayTriangulation(p, constraints);
G      = triangleGrid(delTri.Points, delTri.ConnectivityList);
G      = computeGeometry(G);
%% Add fracture cells into triangle grid 
% Find faces of fracture edges
G.faces.tags = zeros(G.faces.num,1);
faceNodes   = sort(reshape(G.faces.nodes,2,[])',2); %#ok<UDIM>
constraints = sort(delTri.Constraints,2);
for iter = 1 : size(constraints,1)
    fracFace = ismember(faceNodes,constraints(iter,:),'rows');
    G.faces.tags(fracFace) = tags(iter);
end

% Fracture apertures
frac_apers = [1e-4, 1e-4, 5e-4, 5e-4, 1e-3, 1e-3];
    
% Apertures of all faces. Zero aperture for non-fractures
aperture = zeros(G.faces.num,1);
for k = 1:length(frac_apers)
    aperture(G.faces.tags == k) = frac_apers(k);
end

% Add fracture cells into triangle grid
G  = addhybrid(G,G.faces.tags > 0,aperture);

figure;hold on; axis equal off
plotGrid_DFM(G)
plotFractures(G)
drawnow
%% Set rock parameters
% Find indices of hybrid cells
hybridInd = find(G.cells.hybrid);
nCells = G.cells.num;

% Define permeability and porosity
rock.perm = 1 * milli * darcy * ones(nCells,2);
rock.poro = 0.1 * ones(nCells,1);

% Much higher values are used in the fracture
rock.perm(hybridInd,:) = aperture(G.cells.tags(hybridInd)).^2/12 * [1 1];
rock.poro(hybridInd)   = 1;
%% Define wells and simulation schedule
% well coordinates
wcoords = [10, 10; 90, 90];
% well cells
wc = delTri.pointLocation(wcoords);

% add injector
rate = 0.5*meter^3/day;
W = addWell([], G, rock, wc(1), 'Name', 'I1', 'sign',  1,...
    'comp_i', [1 0], 'Val', rate, 'Type', 'rate');

% add producer
bhp = 100*barsa;
W = addWell(W,  G, rock, wc(2), 'Name', 'P1', 'sign', -1,...
    'comp_i', [1 1], 'Val', bhp, 'Type', 'bhp');

% Define the timesteps
timesteps = [ones(20,1)*1*day; ones(20,1)*2*day; ones(50,1)*10*day];

% Set up the schedule containing both the wells and the timesteps
schedule = simpleSchedule(timesteps, 'W', W);
%% Set up simulation model
fluid = initSimpleADIFluid('mu',    [1, 5, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0]);
                       
% Constant oil compressibility
c        = 0.001/barsa;
p_ref    = 300*barsa;
fluid.bO = @(p) exp((p - p_ref)*c);

gravity reset off
model = TwoPhaseOilWaterModel(G, rock, fluid);
% set up dfm operators
model.operators    = operators_dfm(G, rock);
model.outputFluxes = false;
%% Define initial state
state0   = initResSol(G,  300*barsa, [0, 1]);
%% Simulate base case
[wellSols, states, report] = simulateScheduleAD(state0, model, schedule);
%% Plot oil saturation and pressure maps
figure
for k = 1:5:numel(states)
    RT = report.ReservoirTime(k)/day;
    
    subplot(1,2,1)
    cla
    plotCellData_DFM(G, states{k}.s(:,2))
    caxis([0 1])
    axis equal off
    title(sprintf('Oil Saturation,  %d day', RT))
    pause(0.1)
    
    subplot(1,2,2)
    cla
    plotCellData_DFM(G, states{k}.pressure/barsa)
    axis equal off
    title(sprintf('Pressure,  %d day', RT))
    pause(0.1)
    
end
%% Plot injection pressure and prodution rate
%  injection pressure
bhp  = cellfun(@(x)horzcat(x.bhp)/barsa,  wellSols,  'UniformOutput', false);
bhp  = cell2mat(bhp);

%  prodution rate
qWs  = cellfun(@(x)horzcat(x.qWs)*day,  wellSols,  'UniformOutput', false);
qWs  = cell2mat(qWs);

qOs  = cellfun(@(x)horzcat(x.qOs)*day,  wellSols,  'UniformOutput', false);
qOs  = cell2mat(qOs);


RT   = report.ReservoirTime/day;

figure
plot(RT,  bhp(:,1), 'g.-')
grid on
xlim([0 600])
xlabel('Simulation Time (day)')
ylabel('Injetion Pressure (barsa)')
title('Injetion Pressure of I1')


figure; hold on
plot(RT,  -qWs(:,2), 'b.-')
plot(RT,  -qOs(:,2), 'r.-')
grid on
xlim([0 600])
xlabel('Simulation Time (day)')
ylabel('Injetion Pressure (barsa)')
title('Production Rate of P1')
legend({'Water Rate', 'Oil Rate'}, 'Location','NorthWest')

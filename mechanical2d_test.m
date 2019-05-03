    %TEST FLUID2D_SOLVER
%   Part of MILAMIN: MATLAB-based FEM solver for large problems 
%   Version 1.0.1
%   Copyright (C) 2011, M. Dabrowski, M. Krotkiewski, D.W. Schmid
%   University of Oslo, Physics of Geological Processes
%   http://milamin.org
%   See License file for terms of use.

%==========================================================================
% CLEARING AND INITIALIZATION 
%==========================================================================

%CLEAR ENVIRONMENT, BUT NOT BREAKPOINTS
clc; clear variables;

%SET THE DEFAULT ROOT RENDERER TO ZBUFFER,
set(0, 'DefaultFigureRenderer', 'zbuffer');

%==========================================================================
% USER INPUT: PHYSIX
%==========================================================================
D           = [1e-3;  1];                %Viscosity
Rho         = [   1;  2];                %Density
G           = [   0;  0];                %Gravity

%==========================================================================
% MESH GENERATION: 
% no_pts            - number of points on inclusion outline
% radius            - size of inclusion(s)
% (type = 1)        - square box + circular inclusion
% (type = 2)        - square box + circular hole + circular inclusion
% (mode = 'ascii')  - Triangle output written as ASCII files
% (mode = 'binary') - Triangle output written as binary files
%==========================================================================
fprintf(1, 'PREPROCESSING:      '); tic
no_pts =      60;
radius =     0.2;
type   =       2;
mode   = 'ascii'; 
[GCOORD, ELEM2NODE, Point_id, Phases] = generate_mesh(no_pts,radius,type,mode);

nnod    = size(GCOORD,2);
nel     = size(ELEM2NODE,2);

%add 7th node
ELEM2NODE(7,:)  = nnod+1:nnod+nel;
GCOORD          = [GCOORD, [...
    mean(reshape(GCOORD(1, ELEM2NODE(1:3,:)), 3, nel));...
    mean(reshape(GCOORD(2, ELEM2NODE(1:3,:)), 3, nel))]];

nnod    = size(GCOORD,2);

%==========================================================================
% BOUNDARY CONDITION: PURE SHEAR
%==========================================================================
Bc_ind  = find(Point_id==1);
Bc_val  = [GCOORD(1,Bc_ind)  -GCOORD(2,Bc_ind)];
Bc_ind  = [2*(Bc_ind-1)+1       2*(Bc_ind-1)+2];
fprintf(1, [num2str(toc,'%8.6f')]);
fprintf(1, ['\n Number of nodes:   ', num2str(nnod)]);
fprintf(1, ['\n Number of elems:   ', num2str(nel),'\n']);


%==========================================================================
% SOLVER
%  nip                - Number of integration points (6 or higher)
% (reorder = 'amd')   - AMD reordering
% (reorder = 'metis') - METIS reordering
% (method = 'std')    - Standard matrix computation
% (method = 'opt')    - Optimized matrix computation
%==========================================================================
nip      =       6;  
reorder  =   'amd';
method   =   'opt'; 

[V PRESSURE] = mechanical2d(ELEM2NODE, Phases, GCOORD, D, Rho, G, Bc_ind, Bc_val, nip, reorder, method);

%==========================================================================
% POSTPROCESSING
%==========================================================================
fprintf(1, 'POSTPROCESSING:     '); tic
hh = figure(1);
trisurf(reshape(1:3*nel,3, nel)', GCOORD(1,ELEM2NODE(1:3,:)), GCOORD(2,ELEM2NODE(1:3,:)), zeros(size(GCOORD(1,ELEM2NODE(1:3,:)))), PRESSURE/D(1));
title('Pressure'); view(2); shading interp; colorbar; axis image; axis off
fprintf(1, [num2str(toc,'%8.6f'),'\n']);
fprintf(1, ['\n']);

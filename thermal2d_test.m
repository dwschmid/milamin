%TEST THERMAL2D 
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

%SET THE DEFAULT ROOT RENDERER TO ZBUFFER
%(Nicer plots than OpenGL, no transparency support)
set(0, 'DefaultFigureRenderer', 'zbuffer');

%==========================================================================
% PHYSIX
%==========================================================================
D           = [1; 1000];                %Diffusivities

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

%==========================================================================
% BOUNDARY CONDITION: LINEAR TEMPERATURE PROFILE ON BOX, ZERO FLUX ON HOLE
%==========================================================================
Bc_ind  = find(Point_id==1);
Bc_val  = -GCOORD(1,Bc_ind);
fprintf(1, [num2str(toc,'%8.6f')]);
fprintf(1, ['\n Number of nodes:   ', num2str(nnod)]);
fprintf(1, ['\n Number of elems:   ', num2str(nel),'\n']);

%==========================================================================
% SOLVER
%  nip                - Number of integration points (3 or higher)
% (reorder = 'amd')   - AMD reordering
% (reorder = 'metis') - METIS reordering
% (method = 'std')    - Standard matrix computation
% (method = 'opt')    - Optimized matrix computation
%==========================================================================
nip       =         3; 
reorder   =     'amd';
method    =     'opt';

Temp = thermal2d(ELEM2NODE, Phases, GCOORD, D, Bc_ind, Bc_val, nip, reorder, method);

%==========================================================================
% POSTPROCESSING
%==========================================================================
fprintf(1, 'POSTPROCESSING:     '); tic
figure(1);
trisurf(ELEM2NODE(1:3,:)', GCOORD(1,:), GCOORD(2,:), zeros(size(GCOORD(1,:))),Temp,'EdgeColor','k', 'FaceColor', 'interp');
title('Temperature');view(2); colorbar; axis image, axis off
fprintf(1, [num2str(toc,'%8.6f'),'\n']);
fprintf(1, ['\n']);

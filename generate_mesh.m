function [GCOORD,ELEM2NODE, Point_id, Phase_id] = generate_mesh(no_pts_incl, radius, type, mode)

%   Part of MILAMIN: MATLAB-based FEM solver for large problems 
%   Version 1.0.1
%   Copyright (C) 2011, M. Dabrowski, M. Krotkiewski, D.W. Schmid
%   University of Oslo, Physics of Geological Processes
%   http://milamin.org
%   See License file for terms of use.

modelname = 'model';
pointlist   = [];
segmentlist = [];
regionlist  = [];
holelist    = [];

x_min		= -1;
x_max		=  1;
y_min		= -1;
y_max		=  1;

pts_l = 1;
pts_u = 0;

%BOX
BOX          = [x_min x_max x_max x_min; y_min y_min y_max y_max];
no_pts       = size(BOX,2);
pts_u        = pts_u + no_pts;
BOX_s        = [pts_l:pts_u;pts_l+1:pts_u+1; 1 1 1 1];
BOX_s(2,end)   = pts_l;
pts_l        = pts_l+no_pts;
BOX_p        = [x_min+1e-5; y_min+1e-5;1];

pointlist    = [pointlist   BOX];
segmentlist  = [segmentlist BOX_s];
regionlist   = [regionlist  BOX_p];

%CREATE CIRCLE
theta        = linspace(0,2*pi,no_pts_incl);
theta(end)   = [];
xx           = cos(theta);
yy           = sin(theta);

switch type
    case 1
        %INCLUSION
        alpha        = 0.5;
        center_x     = alpha*x_max+(1-alpha)*x_min;
        center_y     = 0.5*(y_max+y_min);
        INCLUSION    = [center_x + radius*xx; center_y + radius*yy];
        no_pts       = size(INCLUSION,2);
        pts_u        = pts_u + no_pts;
        INCLUSION_s        = [pts_l:pts_u;pts_l+1:pts_u+1; 100*ones(1,no_pts)];
        INCLUSION_s(2,end) = pts_l;
        INCLUSION_p  = [center_x; center_y; 2];
        pointlist    = [pointlist   INCLUSION];
        segmentlist  = [segmentlist INCLUSION_s];
        regionlist   = [regionlist  INCLUSION_p];

    case 2
        %INCLUSION
        alpha        = 0.75;
        center_x     = alpha*x_max+(1-alpha)*x_min;
        center_y     = 0.5*(y_max+y_min);
        INCLUSION    = [center_x + radius*xx; center_y + radius*yy];
        no_pts       = size(INCLUSION,2);
        pts_u        = pts_u + no_pts;
        INCLUSION_s        = [pts_l:pts_u;pts_l+1:pts_u+1; 100*ones(1,no_pts)];
        INCLUSION_s(2,end) = pts_l;
        pts_l        = pts_l+no_pts;
        INCLUSION_p  = [center_x; center_y; 2];
        pointlist    = [pointlist   INCLUSION];
        segmentlist  = [segmentlist INCLUSION_s];
        regionlist   = [regionlist  INCLUSION_p];

        %HOLE
        alpha = 0.25;
        center_x     = alpha*x_max+(1-alpha)*x_min;
        center_y     = 0.5*(y_max+y_min);
        HOLE         = [center_x + radius*xx; center_y + radius*yy];
        no_pts       = size(HOLE,2);
        pts_u        = pts_u + no_pts;
        HOLE_s       = [pts_l:pts_u;pts_l+1:pts_u+1; 100*ones(1,no_pts)];
        HOLE_s(2,end)  = pts_l;
        HOLE_p       = [center_x; center_y];
        %
        pointlist    = [pointlist   HOLE];
        segmentlist  = [segmentlist HOLE_s];
        holelist     = [holelist  HOLE_p];

end

no_pts       = size(pointlist,2);
no_seg       = size(segmentlist,2);
no_reg       = size(regionlist,2);
no_hol       = size(holelist,2);

pointlist    = [1:no_pts;pointlist];
segmentlist  = [1:no_seg;segmentlist];
regionlist   = [1:no_reg;regionlist];
holelist     = [1:no_hol;holelist];

%write triangle input file
fid     = fopen('model.poly','w');
fprintf(fid,'%d 2 0 0\n', no_pts);
fprintf(fid,'%d %15.12f %15.12f\n', pointlist);
fprintf(fid,'%d 1\n', no_seg);
fprintf(fid,'%d %d %d %d\n', segmentlist);
fprintf(fid,'%d\n',no_hol);
fprintf(fid,'%d %15.12f %15.12f\n', holelist);
fprintf(fid,'%d 0\n', no_reg);
fprintf(fid,'%d %15.12f %15.12f %d\n', regionlist);
fclose(fid);

area_glob = (2*pi*radius/no_pts_incl)^2;

if(strcmp(mode,'binary'))
    binary_flag = '-b';
else
    binary_flag = [];
end

system(['triangle ',binary_flag,' -pQIo2q33Aa',num2str(area_glob,'%12.12f'),' ',modelname,'.poly']);



if(strcmp(mode,'binary'))

    %NODES READING
    fid=fopen([modelname,'.node'], 'r');
    fseek(fid, 0, 1);
    file_size	= ftell(fid);
    fseek(fid, 0, -1);
    dummy	= fread(fid,file_size/8,'double');
    fclose(fid);

    GCOORD		= [dummy(6:4:end)';dummy(7:4:end)'];
    Point_id	= dummy(8:4:end)';

    %ELEMS READING
    fid=fopen([modelname,'.ele'], 'r');
    fseek(fid, 0, 1);
    file_size	= ftell(fid);
    fseek(fid, 0, -1);
    dummy	= fread(fid,file_size/8,'double');
    fclose(fid);

    ELEM2NODE	= [dummy(5:8:end)';dummy(6:8:end)';dummy(7:8:end)';dummy(8:8:end)';dummy(9:8:end)';dummy(10:8:end)'];
    ELEM2NODE	= int32(ELEM2NODE);
    Phase_id		= dummy(11:8:end)';

else

    %NODES READING
    fid =fopen(strcat(modelname,'.node'),'r');
    tmp = fscanf(fid, '%d',4);
    nnod = tmp(1);
    GCOORD = fscanf(fid, '%e', [4, nnod]);
    fclose(fid);

    GCOORD(1,:)   = [];
    Point_id = GCOORD(end,:);
    GCOORD(end,:) = [];

    %ELEMS READING
    fid =fopen(strcat(modelname,'.ele'),'r');
    tmp = fscanf(fid, '%d',3);
    nel = tmp(1);
    ELEM2NODE = fscanf(fid, '%d',[8, nel]);
    fclose(fid);
    Phase_id  = ELEM2NODE(end,:);
    ELEM2NODE(1,:) = [];
    ELEM2NODE(end,:) = [];
    ELEM2NODE	= int32(ELEM2NODE);

end

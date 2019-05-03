function T = thermal2d(ELEM2NODE, Phases, GCOORD, D, Bc_ind, Bc_val, nip, reorder, method)
% THERMAL2D Two dimensional finite element thermal problem solver of MILAMIN

%   Part of MILAMIN: MATLAB-based FEM solver for large problems, 
%   Version 1.0.1
%   Copyright (C) 2011, M. Dabrowski, M. Krotkiewski, D.W. Schmid
%   University of Oslo, Physics of Geological Processes
%   http://milamin.org
%   See License file for terms of use.

%==========================================================================
% MODEL INFO
%==========================================================================
nnod         = size(GCOORD,2);
nnodel       = size(ELEM2NODE,1);
nel          = size(ELEM2NODE,2);

%==========================================================================
% CONSTANTS
% Adjust nelblo for your specific CPU for optimal performance
%==========================================================================
ndim         =   2;
nelblo       = 760;                                                                                              

%==========================================================================
% BLOCKING PARAMETERS (nelblo must be < nel)
%==========================================================================
nelblo       = min(nel, nelblo);
nblo         = ceil(nel/nelblo);

%==========================================================================
% PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
[IP_X, IP_w] = ip_triangle(nip);                   
[N dNdu]     = shp_deriv_triangle(IP_X, nnodel);   

%==========================================================================
% DECLARE VARIABLES (ALLOCATE MEMORY)
%==========================================================================
K_all        = zeros(nnodel*(nnodel+1)/2,nel); 
Rhs          = zeros(nnod,1);

%==========================================================================
% INDICES EXTRACTING LOWER PART
%==========================================================================
indx_l       = tril(ones(nnodel)); indx_l = indx_l(:); indx_l = indx_l==1;

switch method
    %======================================================================
    % STANDARD VERSION
    %======================================================================        
    case 'std'        
        %==================================================================
        % DECLARE VARIABLES (ALLOCATE MEMORY)
        %==================================================================
        K_elem      = zeros(nnodel,nnodel);   
        
        %==================================================================
        % i) ELEMENT LOOP - MATRIX COMPUTATION
        %==================================================================
        fprintf(1, 'MATRIX ASSEMBLY:    '); tic;
        for iel = 1:nel
            %==============================================================
            % ii) FETCH DATA OF ELEMENT
            %==============================================================
            ECOORD_X = GCOORD(:,ELEM2NODE(:,iel));
            ED       = D(Phases(iel));

            %==============================================================
            % iii) INTEGRATION LOOP
            %==============================================================
            K_elem(:) = 0;
            for ip=1:nip
                %==========================================================
                % iv) LOAD SHAPE FUNCTIONS DERIVATIVES FOR INTEGRATION POINT
                %==========================================================
                dNdui       = dNdu{ip};

                %==========================================================
                % v) CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
                %==========================================================
                J           = ECOORD_X*dNdui;
                detJ        = det(J);
                invJ        = inv(J);

                %==========================================================
                % vi) DERIVATIVES wrt GLOBAL COORDINATES
                %==========================================================
                dNdX        = dNdui*invJ;

                %==========================================================
                % vii) NUMERICAL INTEGRATION OF ELEMENT MATRICES
                %==========================================================
                K_elem      = K_elem + IP_w(ip)*detJ*ED*(dNdX*dNdX');
            end

            %==============================================================
            % ix) WRITE DATA INTO GLOBAL STORAGE
            %==============================================================
            K_all(:,iel)    = K_elem(indx_l);
        end
        fprintf(1, [num2str(toc),'\n']);
                
    %======================================================================
    % OPTIMIZED VERSION
    %======================================================================    
    case 'opt'        
        %==================================================================
        % DECLARE VARIABLES (ALLOCATE MEMORY)
        %==================================================================        
        K_block     = zeros(nelblo,nnodel*(nnodel+1)/2);
        invJx       = zeros(nelblo, ndim);
        invJy       = zeros(nelblo, ndim);
        il          = 1;
        iu          = nelblo;

        %==================================================================
        % i) BLOCK LOOP - MATRIX COMPUTATION
        %==================================================================
        fprintf(1, 'MATRIX ASSEMBLY:    '); tic;
        for ib = 1:nblo
            %==============================================================
            % ii) FETCH DATA OF ELEMENTS IN BLOCK
            %==============================================================
            ECOORD_x = reshape( GCOORD(1,ELEM2NODE(:,il:iu)), nnodel, nelblo);
            ECOORD_y = reshape( GCOORD(2,ELEM2NODE(:,il:iu)), nnodel, nelblo);
            ED       = reshape(D(Phases(il:iu)),nelblo,1);

            %==============================================================
            % iii) INTEGRATION LOOP
            %==============================================================
            K_block(:)  = 0;
            for ip=1:nip
                %==========================================================
                % iv) LOAD SHAPE FUNCTIONS DERIVATIVES FOR INTEGRATION POINT
                %==========================================================
                dNdui       = dNdu{ip};

                %==========================================================
                % v) CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
                %==========================================================
                Jx          = ECOORD_x'*dNdui;
                Jy          = ECOORD_y'*dNdui;
                detJ        = Jx(:,1).*Jy(:,2) - Jx(:,2).*Jy(:,1);

                invdetJ     = 1.0./detJ;
                invJx(:,1)  = +Jy(:,2).*invdetJ;
                invJx(:,2)  = -Jy(:,1).*invdetJ;
                invJy(:,1)  = -Jx(:,2).*invdetJ;
                invJy(:,2)  = +Jx(:,1).*invdetJ;

                %==========================================================
                % vi) DERIVATIVES wrt GLOBAL COORDINATES
                %==========================================================
                dNdx        = invJx*dNdui';
                dNdy        = invJy*dNdui';

                %==========================================================
                % vii) NUMERICAL INTEGRATION OF ELEMENT MATRICES 
                %==========================================================
                weight      = IP_w(ip)*detJ.*ED;

                indx = 1;
                for i = 1:nnodel
                    for j = i:nnodel
                        K_block(:,indx)  =   K_block(:,indx) + ...
                            (dNdx(:,i).*dNdx(:,j)+ dNdy(:,i).*dNdy(:,j)).*weight;
                        indx = indx + 1;
                    end
                end
            end
            %==============================================================
            % ix) WRITE DATA INTO GLOBAL STORAGE
            %==============================================================
            K_all(:,il:iu)	= K_block';
            
            %==============================================================
            % READJUST START, END AND SIZE OF BLOCK. REALLOCATE MEMORY
            %==============================================================
            il  = il+nelblo;
            if(ib==nblo-1)
                nelblo 	= nel-iu;
                K_block	= zeros(nelblo, nnodel*(nnodel+1)/2);
                invJx   = zeros(nelblo, ndim);
                invJy   = zeros(nelblo, ndim);
            end
            iu  = iu+nelblo;
        end
        fprintf(1, [num2str(toc),'\n']);                
end

%==========================================================================
% ix) CREATE TRIPLET FORMAT INDICES
%==========================================================================
fprintf(1, 'TRIPLET INDICES:    '); tic
indx_j = repmat(1:nnodel,nnodel,1); indx_i = indx_j';
indx_i = tril(indx_i); indx_i = indx_i(:); indx_i = indx_i(indx_i>0);
indx_j = tril(indx_j); indx_j = indx_j(:); indx_j = indx_j(indx_j>0);

K_i = ELEM2NODE(indx_i,:); K_i = K_i(:); 
K_j = ELEM2NODE(indx_j,:); K_j = K_j(:);

indx       = K_i < K_j;
tmp        = K_j(indx);
K_j(indx)  = K_i(indx);
K_i(indx)  = tmp;
fprintf(1, [num2str(toc),'\n']);

%==========================================================================
% x) CONVERT TRIPLET DATA TO SPARSE MATRIX
%==========================================================================
fprintf(1, 'SPARSIFICATION:     '); tic
K_all  = K_all(:);
if exist(['sparse2.' mexext], 'file') == 3
    K      = sparse2(K_i, K_j, K_all);
else
    K      = sparse(double(K_i), double(K_j), K_all);
end
clear K_i K_j K_all;
fprintf(1, [num2str(toc),'\n']);

%==========================================================================
% BOUNDARY CONDITIONS
%==========================================================================
fprintf(1, 'BDRY CONDITIONS:    '); tic;
T            = zeros(nnod,1);
T(Bc_ind)    = Bc_val;
Free         = 1:nnod;
Free(Bc_ind) = [];
Rhs          = Rhs -  K*T - K'*T;  
K            = K(Free,Free);
fprintf(1, [num2str(toc),'\n']);

%==========================================================================
% REORDERING
%==========================================================================
fprintf(1, 'REORDERING:         '); tic;
switch reorder
    case 'metis'
        perm = metis(K);
    case 'amd'
        perm = amd(K);
    otherwise
        error('Unknown reordering')
end
fprintf(1, [num2str(toc),'\n']);

%==========================================================================
% FACTORIZATION - ideally L = lchol(K, perm)
%==========================================================================
fprintf(1, 'FACTORIZATION:      '); tic;
K = K(perm,perm);
K = tril(K)+triu(K,1)';
L = chol(K, 'lower');
fprintf(1, [num2str(toc,'%8.6f'),'\n']);

%==========================================================================
% SOLVE
%==========================================================================
fprintf(1, 'SOLVE:              '); tic;
if exist(['cs_ltsolve.' mexext], 'file') == 3 && exist(['cs_lsolve.' mexext], 'file') == 3
    T(Free(perm)) = cs_ltsolve(L,cs_lsolve(L,Rhs(Free(perm))));
else
    T(Free(perm)) = L'\(L\(Rhs(Free(perm))));
end
fprintf(1, [num2str(toc,'%8.6f'),'\n']);

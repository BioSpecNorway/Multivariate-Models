function [Global,Block]=mba(col,nbdim,opt,Y);
%mba - main for mba_algorithm with different deflation procedures
%
%    function [Global,Block]=mba(col,nbdim,option,Y);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%%  Achim Kohler                                                          %
%%                                                                        %
%%  (Achim Kohler:                                                        %
%%   Center for Biospectroscopy and Data Modelling                        %
%%   Matforsk, Norwegian Food Research Institute,                         %
%%   Osloveien 1, 1430 Ås, Norway)                                        %
%%                                                                        %
%%  last update 25.07.07                                                  %
%%                                                                        %
%%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  INPUT
%  col: collection of data
%    (example: col(1)=table1; col(2)=table2; col(3)=table3, where each table is a
%    saisir structure and there is a row to row correspondence between the
%    tables, the tables can have different (numbers of) variables)
%
%  nbdim: number of dimensions
%
%  Y one data block, for ANOVA PLSR the design
%
%  option:
%    option=1 CPCA with orthogonal block loadings (deflation on block loadings)
%    option=2 CPCA classical (the case that corresponds to PCA of X)
%    option=3 CPCA with orthogonal block scores (deflation on block scores)
%    option=4 CPCA classical with orthogonalised block loadings 
%    option=5 MBPLS with orthogonal block loadings
%    option=6 MBPLS classical (the case that corresponds to PCA of X)
%    option=7 MBPLS with orthogonal block scores
%    option=8 MBPLS classical with orthogonalised block loadings 
% 
%
%  OUTPUT:
%  Global: 
%  Block: 
%
%  Both 'Global' and 'Block' contain several 'saisir' structures, with
%  an .info field each to explain the content of the respective saisir
%  structure.
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

%% Check input, check Row to row correspondence
if ((opt<9)&(opt>4)) %MBPLS
   [NY MY]=size(Y.d);
   nBlock=size(col,2);
   for i=1:nBlock
       [NX MX]=size(col(i).d);
       if (NX~=NY)
            ERROR('MSGID', 'No row to row correspondace in input data');
       end
   end
elseif ((opt<5)&(opt>0)) %CPCA
   [NY MY]=size(col(1).d);
   nBlock=size(col,2);
   for i=1:nBlock
       [NX MX]=size(col(i).d);
       if (NX~=NY)
            error('No row to row correspondace in input data');
       end
   end
end

if ((opt<5)&(opt>0))
    [Global,Block]=mba_algorithm(col,nbdim,opt);
end

if ((opt<9)&(opt>4))
    [Global,Block]=mba_algorithmForPLS(col,nbdim,opt,Y);
end







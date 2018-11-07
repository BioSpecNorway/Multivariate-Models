function [Global,Block]=mba_algorithm(col,nbdim,opt,Y);
%mpca_algorithm - multiblock analysis with different deflation procedures
%
%    function [Global,Block]=mba_algorithm(col,nbdim,option,Y);
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
%%  last update 09.12.09                                                  %
%%  030909 improvement of accuracy test for NIPLAS mpca                   %
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
%    option=1 CPCA with orthogonal block loadings
%    option=2 CPCA classical (the case that corresponds to PCA of X)
%    option=3 CPCA with orthogonal block scores
%    option=4 CPCA classical with orthogonalised block loadings 
%    option=5 MBPLS with orthogonal block loadings 
%    option=6 MBPLS classical (the case that corresponds to PCA of X)
%    option=7 MBPLS with orthogonal block scores
%    option=8 MBPLS classical with orthogonalised block loadings 
%        
%           %%%%%%%% ATTENTION: different options for MBPLS do not exist yet!!!!!!
%
%  OUTPUT:
%  Global: containing global results as super scores and super loadings
%  Block: containing block results as block scores and loadings
%
%  Both 'Global' and 'Block' contain several 'saisir' structures, with
%  an .info field each to explain the content of the respective saisir
%  structure.
%
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optdeflation=opt;

if (opt<5)
    % mpca
    mpca=1;
    mpls=0;
elseif (5<=opt<=12)
    % mpca
    mpca=0;
    mpls=1;
end
    

if (opt==5)
    optdeflation=1;
elseif (opt==6)
    optdeflation=2;
elseif (opt==7)
    optdeflation=3;
elseif (opt==8)
    optdeflation=4;
end


%% Calculate means
%% Variables and observations are coded according to the tables
varBlockCode=[];
obsBlockCode=[];
nBlock=size(col,2);
for i=1:nBlock
    aux=col(i);
    [aux ave]=center(aux);   % centers every column
    Block(i).average=ave;
    xnorm(i)=sqrt(sum(sum(aux.d.*aux.d)));  % calculates the norm for every table
    Block(i).norm.d=xnorm(i);
    Block(i).norm.info='Norm for every table';
    [n,p]=size(aux.d);
    varBlockCode=[varBlockCode;ones(p,1)*i];%% each variable is labelled according to the table 
    obsBlockCode=[obsBlockCode;ones(n,1)*i];%% each observation is labelled according to the table
end


%% Calculate means (mean centering is done outside this program
if (mpls==1) %MBPLS
    [n K2]=size(Y.d);
    ybar=mean(Y.d);
    YOrig=Y;
end


%% Build the global X
whole.d=[];
id_group.d=[];
aux1=(1:n)';
nBlock=size(col,2);
for i=1:nBlock
    aux=col(i);
    [n1,K(i)]=size(aux.d);
    id_group.d=[id_group.d;aux1];   % used for barycentre presentation
    whole.d=[whole.d aux.d];            
    variable{i}=aux.v;
end
whole.i=col(1).i;% Names of the first table are taken as observation names;
whole.v=[num2str(varBlockCode) char(variable)];
wholeOrig=whole;  % Save the not deflated global X

% Calculate the Norm of individual X-tables for the explained variances
% Normalisation is done outside the program
deb=1;
xend=0;
for i=1:nBlock
    xend=size(col(i).d,2)+xend;    
    XI(i)=selectcol(whole,deb:xend);
    NormStart(i)=sum(sum(XI(i).d.*XI(i).d));
    deb=xend+1;
end

% Calculate the Norm of whole for the explained variances
NormWholeStart=sum(sum(wholeOrig.d.*wholeOrig.d));


% Define some variables for MBPLSR
if (mpls==1) %MBPLS
   [n K1]=size(whole.d);
   [n K2]=size(Y.d);

   R=eye(K1);
   I=R;
   BET=zeros(K1,K2,nbdim);
   %previous=BET(:,:,1);    %previous is never used? removed it from code(Sahar)
   xbar=mean(whole.d);      % used for calculating BET in PLSR
end



for idim=1:nbdim   % loop for components
    
    if (mpca==1) %cpca
        idim
        [TGlobal,PGlobal]=DecomposeCPCA(whole,K);
        Global.eigenvec.d(:,idim)=PGlobal;
        Global.eigenvec.info='normalised eigenvectors of global PCA';
        Global.score.d(:,idim)=TGlobal;
        
        [XI,Xhut,Block]=DeflateCPCA(whole,PGlobal,TGlobal,K,optdeflation,idim,Block);
           
    elseif (mpls==1)  %MBPLS
        
        [TGlobal,PGlobal,WGlobal]=DecomposeMPLS(whole,Y,K);
        Global.eigenvec.d(:,idim)=PGlobal;
        Global.eigenvec.info='normalised loadings of MPLS';
        Global.score.d(:,idim)=TGlobal;
        Global.eigenvecweights.d(:,idim)=WGlobal;
        Global.eigenvecweights.info='normalised loading weights of MPLS';
            
        [XI,Xhut,Y,Block]=DeflateMPLS(whole,Y,TGlobal,K,idim,Block);
        
        % Regression coefficient calculations
        V=Global.eigenvecweights.d/(Global.eigenvec.d'*Global.eigenvecweights.d);
        Q=Y.d'*Global.score.d;
        BET(:,:,idim)=V*Q'; 
                
    end
    
    for i=1:nBlock
        % Calculate the Norm of ind. X-residual-tables for the explained variances
        Norm(i)=sum(sum(Xhut(i).d.*Xhut(i).d));
        
        % Explained variance for every block
        ExVar(i,idim)=100*Norm(i)/NormStart(i);
    end
        
    % Here the deflation step follows (the deflation is made on whole.d
    wholeOld=whole;
   
    % put the residuals together to one table
    whole.d=[];
    for i=1:nBlock
        aux=XI(i);
        whole.d=[whole.d aux.d];
    end
        
    % Explained variance for the whole matrix
    Diff=wholeOld.d-whole.d;
    NormWhole=sum(sum(Diff.*Diff));
    ExVarGlobal(idim)=100*NormWhole/NormWholeStart;
   
end

%% store the explained variance using the normalised partial loadings
for j=1:nbdim
   chaine=['PC' num2str(j) '             '];
   Global.score.v(j,:)=chaine(1:10);
   CompLabels(j,:)=chaine(1:10);
end
Global.score.v=[Global.score.v num2str(0.1*round(10*ExVarGlobal')) char(ones(nbdim,1)*'%')];
Global.score.i=whole.i;
Global.ExVar.d=ExVarGlobal;
Global.ExVar.v=CompLabels;
Global.ExVar.i='Explained variances in every step substr. the single blocks from X';


if (mpca==1)
    % calculate the global scores
    P=Global.eigenvec.d;   % the eigenvectors are normalised above for all options %%%% NOT ANYMORE!!!!!
    T=wholeOrig.d*P;
    Global.score.d=T;
    Global.score.info='global scores are calculated from normalised and orthogonal global eigenvectors'; 
elseif (mpls==1)
    %% to be filled ???
end

%% store the same Global properties for all options
for j=1:nbdim
   chaine=['PC' num2str(j) '             '];
   Global.eigenvec.v(j,:)=chaine(1:4);
   Global.eigenvecweights.v(j,:)=chaine(1:4);
   if (optdeflation==1)
        Global.eigenvecACOM.v(j,:)=chaine(1:4);
   end
end
Global.eigenvec.i=wholeOrig.v;
Global.eigenvecweights.i=wholeOrig.v;
if (optdeflation==1)
    Global.eigenvecACOM.i=wholeOrig.v;
end


%% store the same block properties for all options
for i=1:nBlock
    Block(i).eigenvec.i=col(i).v;
    Block(i).eigenvec.v=Global.eigenvec.v;
    if (optdeflation==3)
        Block(i).eigenvecPI.i=col(i).v;
        Block(i).eigenvecPI.v=Global.eigenvec.v;

        Block(i).eigenvec.i=col(i).v;
        Block(i).eigenvec.v=Global.eigenvec.v;
        
        Block(i).weights.i=col(i).v;
        Block(i).weights.v=Global.eigenvec.v;
        
    end
    Block(i).score.i=col(i).i;            
    for j=1:nbdim
       chaine=['PC' num2str(j) '             '];
       Block(i).score.v(j,:)=chaine(1:10);
       CompLabels(j,:)=chaine(1:10);
    end
    Block(i).score.v=[Block(i).score.v num2str(0.1*round(10*ExVar(i,:)')) char(ones(nbdim,1)*'%')];
    if (optdeflation==3)
        Block(i).scoreTI.v=Block(i).score.v;
        Block(i).scoreTI.i=Block(i).score.i;
    end
    Block(i).ExVar.d=ExVar(i,:);
    Block(i).ExVar.v=CompLabels;
    Block(i).ExVar.i='Explained variances for every single block';
    Block(i).proj.i=col(i).v;
    
    if (optdeflation==1)
        Block(i).scoreR.v=Block(i).score.v;
        Block(i).scoreR.i=Block(i).score.i;
    end
        
end

if (optdeflation==4)
    
    X=wholeOrig.d;
    [T,Sigma,V]=svd(X,0);
 
    deb=1;
    xend=0;

    for i=1:nBlock
        xend=size(col(i).d,2)+xend; 
        PI=[];  
        PI=V(deb:xend,:);
        AI=Sigma*PI';
        [Rb,Sigmab,Qb]=svd(AI,0);
        Block(i).eigenvec.d=Qb;
        Block(i).eigenvec.info='Orthogonal and normalised individual eigenvectors';
        Block(i).score.d=col(i).d*Qb;
        Block(i).score.info='Scores: Projection on orthogonal and normalised individual eigenvectors';
        Block(i).Rotation.d=Rb;
        Block(i).Rotation.info='Rotation matrix for optdeflation 3';
        deb=xend+1;
    end              
end
    
if (mpls==1) %MBPLS
    
        [n K1]=size(whole.d);
        [n K2]=size(Y.d);
        
        for var=1:K2
            %% Regression coefficients
            temp.d=reshape(BET(:,var,:),K1,nbdim);
            temp.i=wholeOrig.v;
            temp.v=num2str((1:nbdim)');
            res{var}.nom=Y.v(var,:);
            res{var}.BETA=temp;
            
            %% 0-Regression coefficients
            temp1.d=ybar(var)*ones(1,nbdim)-xbar*temp.d;
            temp1.i='Constant';
            temp1.v=temp.v;
            res{var}.BETA0=temp1;
            
            %% Predictions 
            temp3.d=wholeOrig.d*temp.d + ones(n,1)*temp1.d;
            temp3.i=Y.i;
            temp3.v=temp.v;
            res{var}.PREDY=temp3;
            
            %% Root mean square error
            delta=Y.d(:,var)*ones(1,nbdim)-temp3.d;
            temp4.d=sqrt(sum(delta.*delta)/n);
            temp4.i='Root mean square error of calibration';
            temp4.v=temp.v;
            
            res{var}.RMSEC=temp4;
            this_y=selectcol(Y,var);
            temp5=cormap(this_y,res{var}.PREDY);
            temp5.d=temp5.d.*temp5.d;
            res{var}.r2=temp5;
        end
        
        for var=1:K2   % variables of Y
            deb=1;
            xend=0;
            for i=1:nBlock
                xend=size(col(i).d,2)+xend;    
                Block(i).Y{var}.BETA=selectrow(res{var}.BETA,deb:xend);
                Block(i).Y{var}.BETA.i=Block(i).Y{var}.BETA.i(:,2:end);
                deb=xend+1;
            end
        end
        Global.PLSRres=res;
  
end
    
end


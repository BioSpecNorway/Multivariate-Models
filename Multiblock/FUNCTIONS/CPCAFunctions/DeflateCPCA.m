function [XI,Xhut,Block]=DeflateCPCA(whole,PGlobal,TGlobal,K,optdeflation,idim,Block);
%Decompose data different procedures: SVD or NIPALS
%
%    function [XI,Xhut,Block]=DeflateCPCA(whole,PGlobal,TGlobal,K,optdeflation,idim,Block);
%
%
%  INPUT
%
%  whole: collection of data
%    (example: whole(1)=table1; whole(2)=table2; whole(3)=table3, where each table is a
%    saisir structure and there is a row to row correspondence between the
%    tables, the tables can have different (numbers of) variables)
%  PGlobal: Global loadings of CPCA
%  TGlobal: Global scores of CPCA
%  K: Vector containing the number of variables in each block
%  optdeflation: 
%    option=1 CPCA with orthogonal block loadings
%    option=2 CPCA classical (the case that corresponds to PCA of X)
%    option=3 CPCA with orthogonal block scores
%    option=4 CPCA classical with orthogonalised block loadings 
%
%
%  OUTPUT:
%
%  TGlobal: Global scores  
%  PGlobal: Global Loadings
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%%  Sahar Hassani                                                         %
%%                                                                        %
%%  Sahar Hassani:                                                        %
%%   Department of Medical Genetics,                                      %
%%   University of Oslo and Oslo University Hospital, Oslo, Norway        %
%%                                                                        %
%%   NORMENT, KG Jebsen Centre for Psychosis Research,                    %
%%           Oslo University Hospital, Oslo, Norway                       %
%%                                                                        %
%%                                                                        %
%%                                                                        %
%%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



nBlock=size(K,2);
deb=1;
xend=0;

for i=1:nBlock   % loop for tables
    xend=K(i)+xend;
    
    % Choose the corresponding table
    XI(i)=selectcol(whole,deb:xend); % whole is changed in every iteration step
    PI=[];
    PI=PGlobal(deb:xend); % optdeflation 1
    TI=[];
    TI=TGlobal; % optdeflation 2
    XIold=XI(i);
    
    if (optdeflation==1) % CPCA with orthogonal block loadings
        NORM=sqrt(PI'*PI); % block loadings (not normalised)
        PI=PI/NORM;
        TI=XI(i).d*PI;   % here we project on the loadings
        XI(i).d=XI(i).d-TI*PI'; % calculate the residuals
        Xhut(i).d=TI*PI';
        % calculate deflation residuals (the part that is subtracted
        % and not represented by TGlobal)
        TIR=TI-TGlobal*NORM;   % to be elaborated ???
    elseif (optdeflation==2 || optdeflation==4) % CPCA classical
        PI=XI(i).d'*TI;   % super scores (normalised) scores
        XI(i).d=XI(i).d-TI*PI'; % calculate the residuals
        Xhut(i).d=TI*PI';
        
        NORM=sqrt(PI'*PI);
        PI=PI/NORM; % normalise PI
    elseif (optdeflation==3) % CPCA with orthogonal block scores
        NORM=sqrt(PI'*PI);
        PI=PI/NORM;
        
        TI=XI(i).d*PI;   % here we project on the loadings
        TII=XI(i).d*PI; % not normalised block scores
        NORM=sqrt(TI'*TI);
        TI=TI/NORM;
        
        QI=XI(i).d'*TI;   % here we project on the scores
        XI(i).d=XI(i).d-TI*QI'; % deflation of the block X
        
        Xhut(i).d=TI*QI';
        
        NORM=sqrt(QI'*QI);
        QI=QI/NORM;
    end
    
    
    %% store the results for each table
    if (optdeflation==1) % CPCA with orthogonal block loadings
        Block(i).eigenvec.d(:,idim)=PI;   % normalised
        Block(i).eigenvec.info='Normalised orthogonal block loadings (CPCA with orthogonal block loadings)';
        Block(i).score.d(:,idim)=TI;
        % store deflation residuals (the part that is subtracted
        % and not represented by TGlobal)
        Block(i).scoreR.d(:,idim)=TIR; % 
        Block(i).score.info='Non orthogonal block scores (CPCA with orthogonal block loadings)';
    elseif (optdeflation==2)   % CPCA classical
        Block(i).eigenvec.d(:,idim)=PI;   % normalised
        Block(i).eigenvec.info='Normalised non orthogonal block loadings (CPCA classical)'
        TI=XIold.d*PI;
        Block(i).score.d(:,idim)=TI;
        Block(i).score.info='Non orthogonal block scores (CPCA classical)'
    elseif (optdeflation==3)   % CPCA with orthogonal block scores
        
        Block(i).eigenvecPI.d(:,idim)=PI; % orthogonal 1. step (CPCA) block loadings, non reconstructive
        Block(i).eigenvecPI.info='Normalised orthogonal 1. CPCA step block loadings (CPCA with orthogonal block scores)';
        Block(i).eigenvec.d(:,idim)=QI;
        Block(i).eigenvec.info='Normalised non orthogonal 2. CPCA step block loadings (CPCA with orthogonal block scores)';
        
        Block(i).scoreTI.d(:,idim)=TI;
        Block(i).scoreTI.info='Normalised orthogonal block scores based on PI (CPCA with orthogonal block scores)';
        Block(i).score.d(:,idim)=TII;
        Block(i).score.info='Non normalised orthogonal block scores (CPCA with orthogonal block scores)';
        % reconstitution of weights of original variables
        if (idim==1)
            Block(i).weights.d(:,1)=Block(i).eigenvecPI.d(:,1);
        else
            Block(i).weights.d(:,1)=Block(i).eigenvecPI.d(:,1);
            pp=eye(size(Block(i).eigenvecPI.d,1));
            XXi=XI(i).d;
            for ii=2:idim
                wi=Block(i).eigenvecPI.d(:,ii-1);
                pi= Block(i).scoreTI.d(:,ii-1)'*XXi/norm(Block(i).score.d(:,ii-1));
                pp=pp*(eye(size(Block(i).eigenvecPI.d,1))-wi*pi);
                Block(i).weights.d(:,ii)=pp*Block(i).eigenvecPI.d(:,ii);
                XXi=XXi-Block(i).scoreTI.d(:,ii-1)*Block(i).scoreTI.d(:,ii-1)'*XXi;
            end
        end
        Block(i).weights.info='Non normalised weights of the original variables(CPCA with orthogonal block scores)';
        
    end

    deb=xend+1;    
    
end
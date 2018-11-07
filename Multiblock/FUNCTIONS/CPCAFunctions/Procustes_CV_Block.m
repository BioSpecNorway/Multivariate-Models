function [ThatAll,ThatBlock,PhatMinusM]= Procustes_CV_Block(Collection,PAll,Aopt,BlockNum,PersonGroup)
%Procustes Calculations
%returns two different E matrices: E0:0 PC, E: 1:A PCS
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
%%   Dataset_num not used???                                              %
%%                                                                        %
%%                                                                        %
%%                                                                        %
%%                                                                        %
%%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NBlocks=size(Collection,2);
NGroups=size(PersonGroup.g.d,1);


for i=1:NGroups 
    index=find(PersonGroup.d==i);
    
    % delete respective segment
    for j=1:NBlocks
        DeleteCollection(j)=deleterow(Collection(j),index);
    end
        
    % Mean centering for cross-validation ??? put it in separate routine
    for jblock=1:NBlocks
        aux=DeleteCollection(jblock);
        [aux ave]=center(aux);
        XmeanMinusM(jblock).d=ave.d; % this will be used for prediction
        DeleteCollection_centered(jblock)=aux;
    end
    
    [MinusSegmentGlobal,MinusSegmentLocal]=mba(DeleteCollection_centered,Aopt,2,[]);
    
    PMinusM(i)=MinusSegmentGlobal.eigenvec;
    
    %----Normalizing the whole matrix with respect to XminusM
    Xall=[];
    Xall=Collection; % 'raw' collection
    nBlock=size(Xall,2);
    %Eps=eps*10000;
    
    %% weighting loop
    start_index=0;
    for j=1:nBlock
        % Mean centering
        aux=Xall(j);
        [aux]=center_sahar(aux,XmeanMinusM(j).d);   % centers every column ???
        % weighting blocks
        xnorm(i)=sqrt(sum(sum(aux.d.*aux.d)));  % calculates the norm for every table
        aux.d=aux.d/xnorm(i); % divides every table by its norm
        Xall(j).d=aux.d;
    end
    
    Xblock=Xall;    
    Xall=Concad_Blocks(Xall,BlockNum);
    
    
    %-------Procustes Calculations
    H=PMinusM(i).d(:,1:Aopt)'*PAll(:,1:Aopt);
    [U,S,V]=svd(H);
    Rhat(i).d=U*V';  % Rotation (orthogonal Matrix)
    % Rotate P
    PhatMinusM(i).d=PMinusM(i).d(:,1:Aopt)*Rhat(i).d;
    %PhatMinusM(i).d=PMinusM(i).d;
    PhatMinusM(i).v=PMinusM(i).v(1:Aopt,:);
    PhatMinusM(i).i=Xall.v;
    % Rotate T
    % ThatMinusM(i).d=TMinusM(i).d*Rhat(i).d;
    
    % Xall is above centered with respect to the mean of Xminus (DeleteCollection)
    ThatAll(i).d=Xall.d*PhatMinusM(i).d;
    ThatAll(i).i=Xall.i;
    ThatAll(i).v=PhatMinusM(i).v;
    XhatAll(i).d=ThatAll(i).d*PhatMinusM(i).d';
    XhatAll(i).i=Xall.i;
    XhatAll(i).v=Xall.v;
    Xhat.d(i,:)=XhatAll(i).d(i,:);
    
    %% change for block stability plots
    deb=1;
    xend=0;
    for a=1:size(PMinusM(i).d,2)
        deb=1;
        xend=0;
        for j=1:nBlock
            aux=[];
            aux=Xblock(j);
            K=size(aux.d,2);
            xend=K+xend;
            %% calculation of block scores
            PI=[];
            PI=PhatMinusM(i).d(deb:xend,a);
            NORM=sqrt(PI'*PI); % block loadings (not normalised)
            PI=PI/NORM;
            deb=xend+1;
            ThatBlock(i).d(:,a,j)=aux.d*PI;%tb for
        end
        X_Concat=Concad_Blocks(Xblock,BlockNum);
        %% deflation
        TI=X_Concat.d*PhatMinusM(i).d(:,a);
        X_Concat.d=X_Concat.d-TI*PhatMinusM(i).d(:,a)';
        deb=1;
        xend=0;
        for j=1:nBlock                    
            K=size(Collection(j).d,2);
            xend=K+xend;
            Xblock(j).d=X_Concat.d(:,deb:xend);
            deb=xend+1;
        end
    end
    ThatBlock(i).i=Xall.i;
    ThatBlock(i).v=PhatMinusM(i).v; 
    
end









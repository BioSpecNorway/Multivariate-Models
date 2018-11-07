function [E0,E]= Ecal_CV2(Collection,A,Group)
%Ecal_CV2 - calculate residual matrices based on Cross Validation
%
%    returns two different E matrices: E0:0 PC, E: 1:A PCS
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
%%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iCPCA=2;

% Number of blocks and samples
NBlocks=size(Collection,2);
[N]=size(Collection(NBlocks).d,1);
K=0;
for iBlocks=1:NBlocks
    K=K+size(Collection(iBlocks).d,2);
end


% Allocate matrices for Residuals
E=zeros(N,K,A);
E0=zeros(N,K);
NGroups=size(Group.g.d,1);


for iSeg=1:NGroups 
    iSeg
    SegIndex=find(Group.d==iSeg);
    mSeg=size(SegIndex,1);
    
    % Delete respective segment
    for jblock=1:NBlocks
        DeleteCollection(jblock)=deleterow(Collection(jblock),SegIndex);
    end
    
    % Mean centering for cross-validation
    for jblock=1:NBlocks
        aux=DeleteCollection(jblock);
        [aux ave]=center(aux); 
        XmeanMinusM(jblock).d=ave.d; % this will be used for prediction
        DeleteCollection_centered(jblock)=aux;     
    end
    
    % CPCA in cross-validation loop
    [MinusSegmentGlobal,MinusSegmentLocal]=mba(DeleteCollection_centered,A,iCPCA); 
    PMinusM=MinusSegmentGlobal.eigenvec;  % select the loadings
    
    % Select segment m 
    for jblock=1:NBlocks
        XmBlock(jblock)=selectrow(Collection(jblock),SegIndex); % 
    end   
    
    % Mean centering with respect to XmeanMinusM
    IndexB=0;
    for jblock=1:NBlocks
        % Mean centering
        aux=XmBlock(jblock);
        [aux]=center_sahar(aux,XmeanMinusM(jblock).d);   % centers every column ???
        XmBlock(jblock).d=aux.d;
        KB=size(XmBlock(jblock).d,2);
        E0(SegIndex,(IndexB+1):(IndexB+KB))=XmBlock(jblock).d;
        IndexB=IndexB+size(XmBlock(jblock).d,2);
    end
    
    Xm=Concad_Blocks(XmBlock,NBlocks);

    % Not necessary to do this for the whole matrix
    Tmhat.d=Xm.d*PMinusM.d;
    
    for a=1:A
        Xm_hat.d=Tmhat.d(:,1:a)*PMinusM.d(:,1:a)';
        E(SegIndex,:,a)=Xm.d-Xm_hat.d;            
    end
    
end



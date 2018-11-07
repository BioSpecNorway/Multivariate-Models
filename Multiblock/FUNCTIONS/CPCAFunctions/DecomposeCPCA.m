function [TGlobal,PGlobal,S1]=DecomposeCPCA(whole,K);
%Decompose data different procedures: SVD or NIPALS
%
%    function [TGlobal,PGlobal]=DecomposeCPCA(whole);
%
%
%  INPUT
%
%  whole: collection of data
%    (example: whole(1)=table1; whole(2)=table2; whole(3)=table3, where each table is a
%    saisir structure and there is a row to row correspondence between the
%    tables, the tables can have different (numbers of) variables)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SVD=0;
NIPALS=1;
eps=0.000000001;  % for iteration loops
nBlock=size(K,2);

if (SVD)
    % SVD
    [U,S,V]=svd(whole.d,0);
    PGlobal=V(:,1);
    TGlobal=U(:,1);
    S1=S(1,1);
    %Global.eigenval.d(:,idim)=S(1,1);
    %Global.eigenval.info='eigenvalues of PCA(X) for every deflation step';
end

if (NIPALS)
    %%%% CPCA algorithm (NIPALS for CPCA) %%%%
    % start with an arbitrary tT
    tT=mean(whole.d(:,:),2);  % start value
    x=1.0;
    tTold=0;
    pTold=0;
    iNIPALS=0;
    
    while (x>eps)
        T=[];
        deb=1;
        xend=0;
        for i=1:nBlock
            % block analysis
            xend=K(i)+xend;
            XI(i)=selectcol(whole,deb:xend); % Choose the corresponding table
            pb=XI(i).d'*tT/(tT'*tT);
            normp=sqrt(pb'*pb);
            pb=pb/normp;
            tb=XI(i).d*pb;
            
            T=[T tb];
            deb=xend+1;
        end
        % global analysis
        wT=T'*tT/(tT'*tT);
        normw=sqrt(wT'*wT);
        wT=wT/normw;
        tT=T*wT;
        
        % Normalization of T and P
        tT=tT/(sqrt(tT'*tT));
        pT=whole.d'*tT;
        pT=pT/(sqrt(pT'*pT));
        
        % test accuracy (on the normalized) T and P
        iNIPALS=iNIPALS+1;
        if (iNIPALS>1)
            xT=(tT-tTold)'*(tT-tTold);
            xP=(pT-pTold)'*(pT-pTold);
            x=max(xT,xP);
        end
        
        tTold=tT;
        pTold=pT;
        
    end
    %%%% end NIPALS %%%%
    PGlobal=pT;
    TGlobal=tT;
end

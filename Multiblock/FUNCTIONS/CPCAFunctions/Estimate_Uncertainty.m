function [Pval]= Estimate_Uncertainty(Aopt,PhatMinusM,PAll,PersonGroup)
%calculates uncertainty for given coefficient
%the program returns pvalue based on a t-test for each variable
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

M=size(PhatMinusM,2);
NGroups=size(PersonGroup.g.d,1);
S2=zeros(size(PAll,1),Aopt);

for a=1:Aopt    
    for m=1:M
        S2(:,a)=S2(:,a)+(PAll(:,a)-PhatMinusM(m).d(:,a)).^2;
    end
end

K=size(PAll,1);
aux=(K/(K-Aopt-1))*((M-1)/M);
S2=S2*aux;


for a=1:Aopt
    tvalue(:,a)=PAll(:,a)./sqrt(S2(:,a));
end

Pval.i=PhatMinusM(1).i;
Pval.v='Pvalue for a';
Pval.d=zeros(size(tvalue,1),Aopt);

for a=1:Aopt
    for i=1:size(tvalue,1)
        if(tvalue(i,a)>0)
            Pval.d(i,a)=2*(1-tp(tvalue(i,a),NGroups));
        elseif(tvalue(i,a)<0 || tvalue(i,a)==0)
            Pval.d(i,a)=2*tp(tvalue(i,a),NGroups);
        else
            Pval.d(i,a)=1;
        end
    end
end

%Pval.d(size(tvalue,1)+1,:)=1;
%Pval.d(:,Aopt+1)=1;
% figure;
% pcolor(-log10(Pval.d)');
% set(gcf,'Color',[1 1 1]);
% colorbar;

LogRes=Pval;
for a=1:Aopt
    for i=1:size(Pval.d,1)
        if Pval.d(i,a)==0.0
            LogRes.d(i,a)=5;
        else LogRes.d(i,a)=-log10(Pval.d(i,a));
            if LogRes.d(i,a)>5
                LogRes.d(i,a)=5;
            end
        end
    end
end
   


figure;
set(gcf,'Color',[1 1 1]);
plot(LogRes.d(:,1),'color','blue');
axis tight;
xlabel('\bf Variables');
ylabel('\bf -log10(p-value)');






function [ExplainedVariance_Percent]= PC_optimalASD(E0,E,A,option,N0,K0,Collection)
%center				- calculates number of PCs using different criterion
%returns optimal number of PCs
%
%     E0:
%     E:
%     A:
%     option
%     N0
%     K0
%
%
%
%


%%%%%%%%%%%% Aopt=6; % just give some number untile decide about the criterion for Aopt
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

iCPCA=2;

% Determine the size of the multiblock set
[N,K]=size(E0);

PRESS=[];

nBlock=size(Collection,2);

%% Calculate PRESS0 for all blocks separately
start=1;
for b=1:nBlock
    KB=size(Collection(b).d,2);
    PRESS0(b)=sum(sum(E0(:,start:start+KB-1).^2));
    start=start+KB;
end

%% Calculate PRESS0 for all blocks together
PRESS0(nBlock+1)=sum(sum(E0(:,:).^2)); 


%% Calculate PRESS
for a=1:A
    start=1;
    for b=1:nBlock % single blocks
        KB=size(Collection(b).d,2);
        PRESS(b,a)=sum(sum(E(:,start:start+KB-1,a).^2)); 
        start=start+KB;
    end
    PRESS(nBlock+1,a)=sum(sum(E(:,:,a).^2));  % total PRESS for all blocks together
end

if (1)
    %% Calculate C, Partial leverage   %% this part needs revision
    C=zeros(nBlock,A);
    for i=1:(nBlock)
        aux=Collection(i);
        [aux ave]=center(aux);   % auxdummy is not used
        Storeaux(i).d=aux;
        Xmean(i)=ave;
        Collection_centered(i)=aux;
    end
    [Global,Block]=mba(Collection_centered,A,iCPCA);
    P=Global.eigenvec.d;
    for a=1:A
        start=1;
        for b=1:nBlock
            %%% Check this %%%
            C(b,a)=sum(sum((P(start:start+size(Collection(b).d,2)-1,1:a)).^2));
            start=start+size(Collection(b).d,2);
        end
    end
end

%%
if option==1 %Row-wise CV, Hararld idea
    for a=1:A
        for b=1:nBlock
            denominator=((N(1)-N0)*(size(Collection(b).d,2)-K0-C(b,a))); %
            %??????
            KB=size(Collection(b).d,2);
            %denominator=((N(1)-N0)*(KB-K0));
            if denominator==0
                denominator=1;
            elseif denominator<0
                denominator=(-1)*denominator;
            end                
            MPRESS(b,a)=PRESS(b,a)/denominator; %Mean predicted residual sum of squares
        end
        denominator0=(N(1)-N0)*(K-K0-a);
        if denominator0==0
           denominator0=1; 
        elseif denominator0<0
            denominator0=(-1)*denominator0;            
        end
        MPRESS(nBlock+1,a)=PRESS(nBlock+1,a)/denominator0;
    end
    
    for b=1:nBlock
        denominator=(N(1)-N0)*(size(Collection(b).d,2)-K0);
        if denominator==0
            denominator=1;
        elseif denominator<0
            denominator=(-1)*denominator;
        end
        MPRESS0(b)=PRESS0(b)/((N(1)-N0)*(size(Collection(b).d,2)-K0)); %Mean predicted residual sum of squares
    end
    MPRESS0(nBlock+1)=PRESS0(nBlock+1)/((N(1)-N0)*(K-K0));
    
    
    
    for b=1:nBlock+1 
        MSE(b,1)=MPRESS0(b);
        MSE(b,2:A+1)=MPRESS(b,1:A);
    end
    
    %% correction for Plots, Harald's idea
    
        for a=1:A
            for b=1:nBlock+1
                MSE(b,a)=MSE(b,a)+(a-1)*MSE(b,1)*0.03;
            end
        end
    
    %%
    
    RMSE=sqrt(MSE);
    
    aux=MSE(:,1);    
    ExplainedVariance=(aux*ones(1,size(MSE,2))-MSE)./(aux*ones(1,size(MSE,2)));
    ExplainedVariance2=ExplainedVariance(:,2:A+1).*100;
  
    figure;%RMSE Line plot
    subplot(1,2,1);
    RMSE_max=RMSE((nBlock+1),1);
    RMSE_norm(nBlock+1,:)=RMSE((nBlock+1),1:A+1)/(RMSE_max);
    axis tight;
    col={[0 1 1],[0 0 1],[0 0 0],[1 0 1],[1 0 0],[0.5 0.5 0.5],[0 1 1],[0 0 1],[0 0 0],[1 0 1],[1 0 0],[0.5 0.5 0.5]};
    sx=0.1;
    
    if 0
        hg=plot(0:A-1,RMSE_norm(nBlock+1,1:A),':r');
        hold on;
        plot(0:A-1,RMSE_norm(nBlock+1,1:A),'.r');
        set(gcf,'Color',[1 1 1]);
        for x=0:A-1
            text((x+0.1),RMSE((nBlock+1),(x+1))/(RMSE_max),'g','color','r');
        end
        xlabel('PC number');
        ylabel('Root Mean Square Error');
        for b=1:nBlock
            RMSE_max=max(RMSE(b,1));
            RMSE_norm(b,:)=RMSE(b,1:A+1)/(RMSE_max);
            h(b)=plot(0:A-1,RMSE_norm(b,1:A),'k','color',col{b});
            hold on;
            plot(0:A-1,RMSE_norm(b,1:A),'.','color',col{b});
            for x=0:A-1
                text((x+sx),RMSE(b,(x+1))/(RMSE_max),num2str(b),'color',col{b});
            end
        end
        legend([hg,h],{'Global','TimeBlock','FunctionBlock','ComorbiditiesBlock'});
    end
    ExplainedVariance_Percent=ExplainedVariance.*100;
    
    if 1
        %subplot(1,2,2)
        for b=1:nBlock
            h(b)=plot(0:A-1,ExplainedVariance_Percent(b,1:A)','k','color',col{b});
            hold on;
            plot(0:A-1,ExplainedVariance_Percent(b,1:A)','.','color',col{b});
            for x=0:A-1
                text((x+sx),ExplainedVariance_Percent(b,(x+1)),num2str(b),'color',col{b});
            end
            
        end
        hg=plot(0:A-1,ExplainedVariance_Percent(nBlock+1,1:A)',':r');
        hold on;
        plot(0:A-1,ExplainedVariance_Percent(nBlock+1,1:A)','.r');
        for x=0:A-1
            text((x+sx),ExplainedVariance_Percent(nBlock+1,(x+1)),'g','color','r');
        end
        
        xlabel('PC number');
        ylabel('Cross-validated Explained Variance');
        text(0.1,87,'\bf b');
        %title('\bf a');
        legend([hg,h],{'Global','TimeBlock','FunctionBlock','ComorbiditiesBlock'});
        
    end
    
    
   
    
    %% sahar
    
    step=[];
    for block=1:nBlock+1
        for Component=1:A-1 %%%%%%%%%%MODIFY 
            step(block,Component)=ExplainedVariance_Percent(block,Component+1)-ExplainedVariance_Percent(block,Component);
        end
    end
    
    %figure;
    subplot(1,2,2);
    set(gcf,'Color',[1 1 1]);
    bar(step',1);
    %colormap Lines;
    xlabel('PC number');
    ylabel('Cross-validated Explained Variance by each PC');
    legend('TimeBlock','FunctionBlock','ComorbiditiesBlock','Global');
    
    %%
    
    
elseif option==2%Validation Method that used in CV of Wold
    
    [Collection_centred ave]=center(Collection);
    [Global,Block]=mba_sahar_test(Collection_centred,A,2);
    P=Global.eigenvec.d;
    T=Global.score.d;

    SS0=sum(sum(Collection_centred.d(:,:).^2));
    for a=1:A
        Collection_centred.d=Collection_centred.d-T(:,a)*P(:,a)';
        SS(a+1)=sum(sum(Collection_centred.d(:,:).^2));
    end
    %SS(1)=SS0;
    R=PRESS(1)./SS;
    
    Y_min=min(min(R(:,:)));
    Y_max=max(max(R(:,:)));
    figure;
    h=plot(1:A,R(1:A));
    set(h,'LineWidth',1);
    title(texlabel('R=PRESS./SS'));
    grid on;
    figure;
    Y_min=min(min(R(:,:)));
    Y_max=max(max(R(:,:)));
    h=plot(2:12,R(2:12));
    set(h,'LineWidth',1);
    title(texlabel('R=PRESS./SS'));
    grid on;
elseif option==3
    for b=1:nBlock
        figure;
        Y_min=min(min(PRESS(b,:)));
        Y_max=max(max(PRESS(b,:)));
        h=plot(1:12,PRESS(b,1:12));
        set(h,'LineWidth',1);
        title(texlabel('PRESS(b,a)=sum(sum(E(:,:,a).^2))'));
        grid on;
    end
end
    
    
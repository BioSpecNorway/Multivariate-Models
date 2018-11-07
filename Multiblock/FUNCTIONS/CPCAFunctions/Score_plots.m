function []= Score_plots(Global,Block,PC1,PC2,gr1,gr2)
%				Plottng block score plots & global score plot
%
%
Nb=size(Block,2);

if Nb==6
    figure;
    for i=1:6
        subplot(2,3,i);
        carte_couleur1(Block(i).score,PC1,PC2,gr1,gr2);
        set(gcf,'Color',[1 1 1]);
        FontSize=18;
        title1=['Block ',num2str(i)];
        title(title1);
    end
    
    figure;
    carte_couleur1(Global.score,PC1,PC2,gr1,gr2);
    set(gcf,'Color',[1 1 1]);
    FontSize=18;
    title({'\bf Global'});
    
elseif Nb==5
    figure;
    for i=1:5
        subplot(2,3,i);
        carte_couleur1(Block(i).score,PC1,PC2,gr1,gr2);
        set(gcf,'Color',[1 1 1]);
        FontSize=18;
        title1=['Block ',num2str(i)];
        title(title1);
    end
    
    subplot(2,3,6);
    carte_couleur1(Global.score,PC1,PC2,gr1,gr2);
    set(gcf,'Color',[1 1 1]);
    FontSize=18;
    title({'\bf Global'});
    
elseif Nb==4 
    if 0
    figure;
    for i=1:4
        subplot(2,2,i);
        carte_couleur1(Block(i).score,PC1,PC2,gr1,gr2);
        set(gcf,'Color',[1 1 1]);
        FontSize=18;
        title1=['Block ',num2str(i)];
        title(title1);
    end
    
    figure();
    carte_couleur1(Global.score,PC1,PC2,gr1,gr2);
    set(gcf,'Color',[1 1 1]);
    FontSize=18;
    title({'\bf Global'});
    
    end
    
    
    figure;
    for i=1:2
        subplot(2,3,i);
        carte_couleur1(Block(i).score,PC1,PC2,gr1,gr2);
        set(gcf,'Color',[1 1 1]);
        FontSize=18;
        title1=['Block ',num2str(i)];
        title(title1);
    end
    for i=3:4
        subplot(2,3,i+1);
        carte_couleur1(Block(i).score,PC1,PC2,gr1,gr2);
        set(gcf,'Color',[1 1 1]);
        FontSize=18;
        title1=['Block ',num2str(i)];
        title(title1);
    end
    subplot(2,3,6);    
    carte_couleur1(Global.score,PC1,PC2,gr1,gr2);
    set(gcf,'Color',[1 1 1]);
    FontSize=18;
    title({'\bf Global'});
       
    
    
    
elseif Nb==3
    figure;
    for i=1:3
        subplot(2,2,i);
        carte_couleur1(Block(i).score,PC1,PC2,gr1,gr2);
        set(gcf,'Color',[1 1 1]);
        FontSize=18;
        title1=['Block ',num2str(i)];
        title(title1);
    end
    
    subplot(2,2,4);
    carte_couleur1(Global.score,PC1,PC2,gr1,gr2);
    set(gcf,'Color',[1 1 1]);
    FontSize=18;
    title({'\bf Global'});
elseif Nb==2
    figure;
    for i=1:2
        subplot(1,3,i);
        carte_couleur1(Block(i).score,PC1,PC2,gr1,gr2);
        set(gcf,'Color',[1 1 1]);
        FontSize=18;
        title1=['Block ',num2str(i)];
        title(title1);
    end
    
    subplot(1,3,3);
    carte_couleur1(Global.score,PC1,PC2,gr1,gr2);
    set(gcf,'Color',[1 1 1]);
    FontSize=18;
    title({'\bf Global'});
    
elseif Nb==1
    figure;
    carte_couleur1(Global.score,PC1,PC2,gr1,gr2);
    set(gcf,'Color',[1 1 1]);
    FontSize=18;
    title({'\bf Global'});
end










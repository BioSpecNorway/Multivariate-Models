% Main - main for ASD article CPCA analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%%  Sahar Hassani                                                         %
%%                                                                        %
%%  Sahar Hassani:                                                       %
%%   Department of Medical Genetics,                                      %
%%   University of Oslo and Oslo University Hospital, Oslo, Norway        %
%%                                                                        %
%%   NORMENT, KG Jebsen Centre for Psychosis Research,                    %
%%           Oslo University Hospital, Oslo, Norway                       %
%%                                                                        %
%%  last update 15.10.18                                                  %
%%                                                                        %
%%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Main script
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
addpath(strcat('..\SaisirFunctions'));
addpath(strcat(pwd,'\FUNCTIONS\CPCAFunctions'));

%% Reading the data 

[Data1,Var1]=xlsread('DATA\Data.xlsx',1);
[Data2,Var2]=xlsread('DATA\Data.xlsx',2);
[Data3,Var3]=xlsread('DATA\Data.xlsx',3);


%% Making data blocks in saisir structure
Col(1).d=Data1; %
Col(1).v=char(Var1);
Col(1).i=num2str((1:size(Data1,1))');

Col(2).d=Data2; %
Col(2).v=char(Var2);
Col(2).i=Col(1).i;

Col(3).d=Data3; %
Col(3).v=char(Var3);
Col(3).i=Col(1).i;


% Order Saisir structure
X(1)=ordersaisirstructure(Col(1));
X(2)=ordersaisirstructure(Col(2));
X(3)=ordersaisirstructure(Col(3));

NBlocks=size(X,2);


%% Mean centering & scaling the X data blocks  
% Mean centering
for i=1:(NBlocks)
    aux=X(i);
    [aux ave]=center(aux);   % auxdummy is not used
    Storeaux(i).d=aux;
    Xmean(i)=ave;
    X_centered(i)=aux;   
end

% Scaling data blocks
for i=1:NBlocks
    xnorm(i)=sqrt(sum(sum(X_centered(i).d.*X_centered(i).d)));  % calculates the norm for every table
    X_centered(i).d=X_centered(i).d/xnorm(i); % divides every table by its norm   
end

%% Running CPCA
PC=10;
[Global,Block]=mba(X_centered,PC,2); % calculates 10 PCs. NB! opt=2 is used for running classical CPCA
Score_plots(Global,Block,1,2,1,3); % Score plots for PC1 & PC2 NB! 1,3 are the indeces used for colouring the samples 
Score_plots(Global,Block,3,4,1,3); % Score plots for PC3 & PC4
Score_plots(Global,Block,5,6,1,3); % Score plots for PC5 & PC6

figure;
set(gcf,'Color',[1 1 1]);
correlation_plot_cpca(Global.score,1,2,X_centered(1),X_centered(2),X_centered(3)); % Correlation loading plots for PC1 & PC2

figure;
set(gcf,'Color',[1 1 1]);
correlation_plot_cpca(Global.score,3,4,X_centered(1),X_centered(2),X_centered(3)); % Correlation loading plots for PC3 & PC4

figure;
set(gcf,'Color',[1 1 1]);
correlation_plot_cpca(Global.score,5,6,X_centered(1),X_centered(2),X_centered(3)); % Correlation loading plots for PC5 & PC6


%% Save results
filename = 'RES\GlobalScoresPC1to10.xlsx';
xlswrite(filename,Global.score.d(:,:))

filename = 'RES\BlockScoresPC1to10.xlsx';
xlswrite(filename,Block(1).score.d,1);
xlswrite(filename,Block(2).score.d,2);
xlswrite(filename,Block(3).score.d,3);

filename = 'RES\GlobalLoadingsPC1to10.xlsx';
xlswrite(filename,Global.eigenvec.d(:,:));
xlswrite(filename,Global.eigenvec.i,2);

filename = 'RES\BlockLoadingsPC1to10.xlsx';
xlswrite(filename,Block(1).eigenvec.d,1);
xlswrite(filename,Block(2).eigenvec.d,2);
xlswrite(filename,Block(3).eigenvec.d,3);


%% Cross-validation 
% Calculate residuals based on Cross Validation
Group=create_group1(X_centered(1),1,3);

PC=11;
[E0,E]= Ecal_CV2(X,PC,Group);

% Find the number of components
N0=0; K0=0;
CVExplainedVariance=PC_optimalASD(E0,E,PC-1,1,N0,K0,X);
xlswrite('RES\CVExplainedVariance.xlsx',CVExplainedVariance);









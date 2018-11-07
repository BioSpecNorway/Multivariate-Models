function handle=correlation_plot_cpca(scores,col1, col2, varargin);
%correlation_plot       - Draw a correlation between scores and tables
%function handle=correlation_plot(scores,col1,col2, saisir1,saisir2, ...);
%% Input argument:
%% scores:                  - ORTHOGONAL scores obtaiend by multidimensional analysis 
%% col1 and col2:           - Indices of the scores to be plotted
%% saisir1, saisir2, ...    - Arbitrary number of tables giving the variables to be plotted
%% The number of rows in the scores and other tables must be identical
charsize1=8;
charsize2=8;
charsize3=8;
mSize=8;
couleur=[0 0 0; 0 0 1; 0 0.7 0; 1 0 0; 0.5 0.5 0; 0.5 0 0.5; 0 0.5 0.5 ; 0.25 0.25 0.25 ; 0.5 0 0; 0 0.5 0; 0 0.5 0; 0.1 0.2 0.3; 0.3 0.2 0.1; 0.5 0.5 0.8; 0.1 0.8 0.1 ];  

[n,p]=size(scores.d);

axis square;
t=0:pi/20:2*pi;
plot(sin(t),cos(t),'k')
hold on;
cst=sqrt(0.5);%% drawing 50% variance explained
h=plot(sin(t)*cst,cos(t)*cst,'k');
set(h,'LineStyle','--');
line([-1 1],[0 0],'color','k');
line([0 0],[-1 1],'color','k');
axis square
xlabel(scores.v(col1,:));
ylabel(scores.v(col2,:));
[a,ntable]=size(varargin);
scores=selectcol(scores,[col1 col2]);
scores
cor=cormap(scores,scores);
aux=cor.d(1,2);
if((aux*aux)>0.05)
    disp([ 'Determination coefficient between the 2 scores = ' num2str(aux*aux)]);
    disp('Warning : the scores are not sufficiently orthogonal to allow a correct representation');
end
for i=1:ntable
    this_color=couleur(i,:);
    this_table=(varargin{i});
    cor=cormap(this_table,scores);
    [n,p]=size(cor.d);
    for j=1:n
        Vertnum=(cor.i(j,:));
        text(cor.d(j,1),cor.d(j,2),['  ',(Vertnum)],'FontSize',charsize1,'Color',this_color);
        plot(cor.d(j,1),cor.d(j,2),'.','MarkerSize',mSize,'Color',this_color);
    end    
end

hold off; 

handle=h;

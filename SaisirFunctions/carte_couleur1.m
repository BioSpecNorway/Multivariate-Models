function carte_couleur1(s,col1,col2,startpos,endpos,labcol1,labcol2,titre,charsize,marg)
%carte_couleur1			- colored map : using a portion of the identifiers as labels 
%carte_couleur1(X,col1,col2,startpos,endpos,(col1label),(col2label),(title),(charsize),(margin))
%Biplot of two columns as colored map 
%The coloration of the displayed descriptors depends on the arguments
%startpos and endpos.
%From the names of individual, the string name(sartpos:endpos) is extracted. Two observations
% for which these strings are different, are alos colored differently.   

couleur=[0 0 0; 1 0 0; 0 0 1; 0 0.7 0; 0.5 0.5 0; 0.5 0 0.5; 0 0.5 0.5 ; 0.25 0.25 0.25 ; 0.5 0 0; 0 0.5 0; 0 0.5 0; 0.1 0.2 0.3; 0.3 0.2 0.1; 0.5 0.5 0.8; 0.1 0.8 0.1 ;
         0 0 0.2; 1 0 0.2; 0 0 0.8; 0 0.7 0.2; 0.5 0.5 0.2; 0.5 0 0.7; 0 0.5 0.7 ; 0.25 0.25 0.45 ; 0.5 0 0.2; 0 0.5 0; 0 0.5 0.2; 0.1 0.2 0.2; 0.3 0.2 0.3; 0.5 0.5 1.0; 0.1 0.8 0.3];
cla;


title('');   
margin=0.05;
if(nargin>9) margin=marg;end;

csize=8;
if(nargin>8)csize=charsize;end;
axis('auto');
%axis([min(s.d(:,col1))*1.1 max(s.d(:,col1))*1.1 min(s.d(:,col2))*1.1 max(s.d(:,col2))*1.1]); 
minx=min(s.d(:,col1));maxx=max(s.d(:,col1));miny=min(s.d(:,col2));maxy=max(s.d(:,col2));
minx=minx-(maxx-minx)*margin;
miny=miny-(maxy-miny)*margin;
maxx=maxx+(maxx-minx)*margin;
maxy=maxy+(maxy-miny)*margin;

axis([minx maxx miny maxy]); 

model(1,:)=s.i(1,startpos:endpos);
nmodel=1;

gr=create_group1(s,startpos,endpos);
indice=gr.d;
for i=1:size(s.i,1)
    text(s.d(i,col1),s.d(i,col2),s.i(i,startpos:endpos),'FontSize',csize,'Color',couleur(mod(indice(i),15)+1,:),'FontSize',10);
end   
if(nargin >5) xlabel(labcol1)
	else xlabel(s.v(col1,:),'FontSize',10);   
end;
if(nargin>6) ylabel(labcol2);
	else ylabel(s.v(col2,:),'FontSize',10);   
end;
if(nargin>7) title(titre);end;
%hold off;

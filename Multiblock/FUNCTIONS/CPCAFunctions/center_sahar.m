function [saisir] = center_sahar(saisir1,ave)
%center				- subtracts the average row to each row
% function [X xmean] = center(X1)
saisir.v=saisir1.v;
saisir.i=saisir1.i;
%size(ones(size(saisir1.d,1),1)*aux)
saisir.d=saisir1.d-ones(size(saisir1.d,1),1)*ave;

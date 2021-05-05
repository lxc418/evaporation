%DATASET 17 locate each node for calculating XX & YY
%%the paprameters below can be got from .INP
NN=2196; %total nodes
NN1=61;  %along x-axis 
NN2=36;  %along y-axis 

%% along y
aa=zeros (NN,3);
aa(:,1)= -(1:NN);
aa(:,2)= mod(-aa(:,1),NN2);
[ii,jj]=find(~-aa(:,2));
aa(ii,2)= NN2;
aa(:,2)=NN2+1-aa(:,2);
aa(:,3)=(-aa(:,1)+aa(:,2)-1)./NN2;

%% along x
% aa=zeros (NN,3);
% aa(:,1)= -(1:NN);
% 
% aa(:,3)= mod(-aa(:,1),NN1);
% [ii,jj]=find(~-aa(:,3));
% aa(ii,3)= NN1;
% 
% aa(:,2)=36-(-aa(:,1)-aa(:,3))./NN1;

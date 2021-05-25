%DATASET 17 locate each node for calculating XX & YY
%%the paprameters below can be got from .INP
clear

NN1=61;  %along x-axis 
NN2=61;  %along y-axis 

%% along y
NN=NN1*NN2;
aa=zeros (NN,3);
i=(1:NN2:NN);
j=(NN1:-1:1);
k=0:NN1-1;
% aa(i,1)=j;
for ii=0:NN2-1
aa(i+ii,1)=-(NN1*(ii+1)-k);
end
aa(:,2)= mod(-aa(:,1),NN1);
[ii,jj]=find(~-aa(:,2));
aa(ii,2)= NN1;
aa(:,2)=NN1+1-aa(:,2);
aa(:,3)=(-aa(:,1)+aa(:,2)-1)./NN1;

%% along x
% aa=zeros (NN,3);
% aa(:,1)= -(1:NN);
% 
% aa(:,3)= mod(-aa(:,1),NN2);
% [ii,jj]=find(~-aa(:,3));
% aa(ii,3)= NN2;
% 
% aa(:,2)=36-(-aa(:,1)-aa(:,3))./NN2;

%  Comparision against the different senarios
lw    =2;lw2         =1;tfs           =10;fs  =12;ms   =3;por=0.39;
tma   =[60;80;100];rh=[0.05;0.01;0];mw=0.018;R=8.314;iv=10;
kg  =5e-12; % permeability of air (m2) by thomas
miua=1.846e-5; % viscosity of air (kg/m/s) by thomas
dz  =0.1;    % meter
dz1 =0.3;   %
dz2 =0.1;
dtk =1;
drh =0.02;
det =0.1;  % the difference of et
%% get the elevation of each sensor
% (E)levation of (S)ensors from DT(8)5g at channel (C) (unit m) for moisture content
% (E)levation of (S)ensors from DT(8)5g at channel (1) (unit m) for temperature and pore water pressure
   es8c=-0.01;es8d=-0.05;es8e=-0.01;es8f=-0.05;es81=-0.03;es82=-0.08;es83=-0.03;es84=-0.08;
   % (E)levation of (S)ensors from EM(5)0 at Channel (A) 1-a 2-b 3-c 4-d ... (unit m)
   es51=0.00;es52=0.1;es53=0.01;es54=0.02;es55=0.02;
   % (E)levation at the (B)ottom of the (9)0mm (C)olumn
   eb9c=-0.2;
   % elevation of the light
% ALL THE UNIT SHOULD BE IN STANDARD UNIT (M KG K)
% Change all of unit of the average temperature values from celsius to kelvin
% at5ak=at51+az;at5bk=at52+az;
% at5ck=at53+274.15;at5dk=at54+az;at5ek=at55+az;                  
% %% Interpolate the temperature profile of the column
%  % (F)itting (T)emperature (P)rofile 
%  ftp=zeros(size(wl),2);
%  for i=1:size(wl,1)
%  ftp(i,:)=polyfit([es5c,es5b,es5a],[at5ck(i)+1,at5bk(i),at5ak(i)],1); 
%  end
%   % (F)itting (L)ower (T)emperature (P)rofile 
%    fltp=zeros(size(wl),2);
%   for i=1:size(wl,1)
%  fltp(i,:)=polyfit([es5a,ebc],[at5ak(i),atbk(i)],1); 
%  end
%  
%  
%  
%  % (F)itting (T)emperature (P)rofile across the whole (C)olumn
%  ftpc=zeros(size(wl),3);
%  for i=1:size(wl,1)
%  ftpc(i,:)=polyfit([es5c,es5b,es5a,ebc],[at5ck(i)+1,at5bk(i),at5ak(i),atbk(i)],2); 
%  end
% 
%   z=-0.6:0.01:0.02; % the unit of the z is meters
%   z2=-0.6:0.01:-0.08;
%   z1=-0.08:0.01:0;
%% calculate surface resistance see page 134 of the blue book
% the time series of surface resistance is based on the dt85g time series
% (A)ccumulative (E)vapora(T)ion of (9)0 mm (C)olumn rate with respected to scale data converted to DT(8)5G time data 
aet9c8  =zeros(nr8,1);
% (T)ransient (E)vapora(T)ion of 90 mm column rate with respected to (S)cale data converted to(2) DT(8)5G time data 
et8     =zeros(nr8,1);
% (E)vapora(T)ion rate calculated through the (R)elative (H)umidity sensor at the tunnel using time series of DT(8)5G
etrh8   =zeros(nr8,1);
etrh812 =zeros(nr8,1);
% (R)elative (H)umidity at EM(5)0 Channel (2) interpolated into the time series of DT(8)5G 
rh518   =zeros(nr8,1);
rh528   =zeros(nr8,1);
rh558   =zeros(nr8,1);
% (T)emperature at EM(5)0 Channel (2) interpolated into the time series of DT(8)5G 
t518    =zeros(nr8,1);
t528    =zeros(nr8,1);
t558    =zeros(nr8,1);

va      =zeros(nr8,1);
va1     =zeros(nr8,1);
va2     =zeros(nr8,1);
advm    =zeros(nr8,1);


% rac,rscc, is the changing resistances see page 203
racc   =zeros(nr8,1);
uracc  =zeros(nr8,1);
uraccet=zeros(nr8,1);
rscc   =zeros(nr8,1);
drva   =zeros(nr8,1);
% (RS) using (T)emperature at EM(5)0, channel (1)
ra     =zeros(ns,4);
rs     =zeros(nr8,4);

rab    =zeros(ns,4);
rsb    =zeros(nr8,4);

% (A)erodynamic (R)esistance using (T)emperature at EM(5)0, channel (1)
% (A)erodynamic (R)esistance using (H)alf (S)aturation
rahs    =zeros(ns,1);
rshs    =zeros(nr8,4);
 
  for i=1:ns   % number of sessions
  if (sc(i)==2)  % calculation
       
      
 %  interpolation to DT(8)5G time series
	aet9c8(sst8(i):set8(i))=interp1(dtsd(ssts(i):sets(i)),aet9c(ssts(i):sets(i)),dt8d(sst8(i):set8(i)));
	et8(sst8(i):set8(i))   =interp1(dtsd(ssts(i):sets(i)),etsgf9(ssts(i):sets(i)),dt8d(sst8(i):set8(i)));
	% the first value is not interpolated so one has to add it manually
	et8(sst8(i))           =et8(sst8(i)+1);
	rh518(sst8(i):set8(i)) =interp1(dt5d(sst5(i):set5(i)),rh51(sst5(i):set5(i)),dt8d(sst8(i):set8(i)));  
	rh528(sst8(i):set8(i)) =interp1(dt5d(sst5(i):set5(i)),rh52(sst5(i):set5(i)),dt8d(sst8(i):set8(i)));  
	rh558(sst8(i):set8(i)) =interp1(dt5d(sst5(i):set5(i)),rh55(sst5(i):set5(i)),dt8d(sst8(i):set8(i)));   
	t518(sst8(i):set8(i))  =interp1(dt5d(sst5(i):set5(i)),t51(sst5(i):set5(i)),dt8d(sst8(i):set8(i)));
	t528(sst8(i):set8(i))  =interp1(dt5d(sst5(i):set5(i)),t52(sst5(i):set5(i)),dt8d(sst8(i):set8(i)));
	t558(sst8(i):set8(i))  =interp1(dt5d(sst5(i):set5(i)),t55(sst5(i):set5(i)),dt8d(sst8(i):set8(i)));
	
    
%     h=figure;i=4;
%     plot(dt5d(sst5(i):set5(i))-dt5d(sst5(i)),t52(sst5(i):set5(i)),'r','LineWidth',4);hold on
%     plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),t528(sst8(i):set8(i)));hold on
%     
%     plot(dt5d(sst5(i):set5(i))-dt5d(sst5(i)),t51(sst5(i):set5(i)),'r','LineWidth',4);hold on
%     plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),t518(sst8(i):set8(i)));hold on
    
    
    
%     h=figure;i=4;
%     plot(dt5d(sst5(i):set5(i))-dt5d(sst5(i)),rh52(sst5(i):set5(i)),'r','LineWidth',4);hold on
%     plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),rh528(sst8(i):set8(i)));hold on
%     
%     plot(dt5d(sst5(i):set5(i))-dt5d(sst5(i)),rh51(sst5(i):set5(i)),'r','LineWidth',4);hold on
%     plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),rh518(sst8(i):set8(i)));hold on
    
% check the interpolation
%     h=figure;
%     plot(dt5d(sst5(i):set5(i))-dt5d(sst5(i)),rh52(sst5(i):set5(i)),'-','LineWidth',4);hold on
%     plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),rh528(sst8(i):set8(i)));hold on
%     h=figure;i=2;
%     plot(dtsd(ssts(i):sets(i))-dtsd(ssts(i)),aet9c(ssts(i):sets(i)),'-','LineWidth',4);hold on
%     plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),aet9c8(sst8(i):set8(i)));hold on
    
% %  calculating aerodynamic resistance
%    j=1; % rh=0.05, tma=60
	for j=1:3
%  ra with zero point as the first node, rh=0.05, tma=60; j=2, rh=0.01, tma=80 ; j=1, rh=0, tma=100
   ra(i,j)=(func.svp(t518(zp8at(i))+az)*1/(t518(zp8at(i))+az)-func.svp(tma(j)+az)*rh(j)/(tma(j)+az))*mw/R/rholw./etsp8(zp8at(i));

%  calculate surface resistance see page 134
   rs(zp8at(i):set8(i),j)=(func.svp(    t518(zp8at(i):set8(i))+az     )*1./(  t518(zp8at(i):set8(i))+az   )-func.svp(tma(j)+az)*rh(j)/(tma(j)+az)) ...
      *mw/R/rholw./etsp8(zp8at(i):set8(i))-ra(i,j);   
%  (ra) with zero point at the very (b)eginning, rh=0.05, tma=60  [rab]
   rab(i,j)=(func.svp(t518(sst8(i))+az)*1/(t518(sst8(i))+az)-func.svp(tma(j)+az)*rh(j)/(tma(j)+az))*mw/R/rholw./etsp8(sst8(i));
  % rabm(i,j)=(func.svp(t518(sst8(i))+az)*1/(t518(sst8(i))+az)-func.svp(tma(j)+az)*rh(j)/(tma(j)+az))*mw/R/rholw./etsp8(sst8(i));
   rsb(sst8(i):set8(i),j)=(func.svp(    t518(sst8(i):set8(i))+az     )*1./(  t518(sst8(i):set8(i))+az   )-func.svp(tma(j)+az)*rh(j)/(tma(j)+az)) ...
      *mw/R/rholw./etsp8(sst8(i):set8(i))-rab(i,j);  
	  
%  (R)erodynamic (R)esistance from (H)alf (S)aturation; rh=0.05, tma=60
   rahs(i,j)=(func.svp(t518(hsp8(i))+az)*1/(t518(hsp8(i))+az)-func.svp(tma(j)+az)*rh(j)/(tma(j)+az)) ...
       *mw/R/rholw./etsp8(hsp8(i));
%  (S)urface (R)esistance from (H)alf (S)aturation
   rshs(hsp8(i):set8(i),j)=(func.svp(t518(hsp8(i):set8(i))+az)*1./(t518(hsp8(i):set8(i))+az)-func.svp(tma(j)+az)*rh(j)/(tma(j)+az)) ...
       *mw/R/rholw./etsp8(hsp8(i):set8(i))-rahs(i,j);
end   

 % ra with zero point the the first node. rha and tma are based on the second rht sensor
 
   j=4;
   
%    ra(i,j)=(func.svp(t518(zp8at(i))+az)*1/(t518(zp8at(i))+az)-func.svp(t528(zp8at(i))+az)*rh528(zp8at(i))/(t528(zp8at(i))+az))*mw/R/rholw./etsp8(zp8at(i));
%  %  calculate surface resistance see page 134
%    rs(zp8at(i):set8(i),j)=(func.svp(    t518(zp8at(i):set8(i))+az     )*1./(  t518(zp8at(i):set8(i))+az   )-func.svp(t528(zp8at(i))+az)*rh528(zp8at(i))/(t528(zp8at(i))+az)) ...
%       *mw/R/rholw./etsp8(zp8at(i):set8(i))-ra(i,j);   
%   
   rab(i,j)=(func.svp(t518(sst8(i))+az)*1/(t518(sst8(i))+az)-func.svp(t528(sst8(i))+az)*rh528(sst8(i))/(t528(sst8(i))+az))*mw/R/rholw./etsp8(sst8(i));
%      rab(sst8(i):set8(i),j)=(func.svp(    t518(sst8(i):set8(i))+az     )*1./(  t518(sst8(i):set8(i))+az   )-func.svp(t528(sst8(i):set8(i)   )+az).*rh528(sst8(i))./(t528(sst8(i):set8(i))+az)) ...
%       *mw/R/rholw./etsp8(sst8(i));
%  calculate surface resistance see page 134
     rsb(sst8(i):set8(i),j)=(func.svp(    t518(sst8(i):set8(i))+az     )*1./(  t518(sst8(i):set8(i))+az   )-func.svp(t528(sst8(i):set8(i)   )+az).*rh528(sst8(i))./(t528(sst8(i):set8(i))+az)) ...
      *mw/R/rholw./etsp8(sst8(i):set8(i))-rab(i,j);
%     rsb(sst8(i):set8(i),j)=(func.svp(    t518(sst8(i):set8(i))+az     )*1./(  t518(sst8(i):set8(i))+az   )-func.svp(t528(sst8(i):set8(i)   )+az).*rh528(sst8(i):set8(i))./(t528(sst8(i):set8(i))+az)) ...
%       *mw/R/rholw./etsp8(sst8(i):set8(i))-rab(i,j);



 
% calculate the evaporation rate based on the relative humidity sensor at the tunnel above the soil surface (m/s)
  etrh8(sst8(i):set8(i))=func.dv(  (t558(sst8(i):set8(i))+t528(sst8(i):set8(i))+2*az)/2     ).*...
      (func.rhovs(t558(sst8(i):set8(i))+az).*rh558(sst8(i):set8(i))-func.rhovs(t528(sst8(i):set8(i))+az).*rh528(sst8(i):set8(i)))...
      /0.1/rholw;
   
  etrh812(sst8(i):set8(i))=func.dv(  (t518(sst8(i):set8(i))+t528(sst8(i):set8(i))+2*az)/2     ).*...
      (func.rhovs(t518(sst8(i):set8(i))+az).*rh518(sst8(i):set8(i))-func.rhovs(t528(sst8(i):set8(i))+az).*rh528(sst8(i):set8(i)))...
      /0.1/rholw;
  
%    aetrh812=func.dv(  (t518(sst8(i):set8(i))+t528(sst8(i):set8(i))+2*az)/2     ).*...
%       (func.rhovs(t518(sst8(i):set8(i))+az).*rh518(sst8(i):set8(i)))...
%       /0.1/rholw*ms2mmd;
%   
%      a2etrh812=func.dv(  (t518(sst8(i):set8(i))+t528(sst8(i):set8(i))+2*az)/2     ).*...
%       (func.rhovs(t518(sst8(i):set8(i))+az))...
%       /0.1/rholw*ms2mmd;
  
    

% air flow between the 1 and 5 sensors  see page notebook 203
  va(sst8(i):set8(i))=-kg/miua/por*( func.svp(t518(sst8(i):set8(i))+az).*rh518(sst8(i):set8(i))- func.svp(t528(sst8(i):set8(i))+az).*rh528(sst8(i):set8(i)) )/dz;
  
  va1(sst8(i):set8(i))=-kg/miua/por*( func.svp(t518(sst8(i):set8(i))+az).*rh518(sst8(i):set8(i))- func.svp(tma(1)+az)*rh(1) )/dz1;
  
  va2(sst8(i):set8(i))=-kg/miua/por*( func.svp(t528(sst8(i):set8(i))+az).*rh528(sst8(i):set8(i))- func.svp(tma(1)+az)*rh(1) )/dz2;
%   adv1(sst8(i):set8(i))=-va(sst8(i):set8(i)).*func.rhovs(t518(sst8(i):set8(i))+az).*rh518(sst8(i):set8(i))  /rholw*ms2mmd;
%   adv2(sst8(i):set8(i))=-va(sst8(i):set8(i)).*func.rhovs(t528(sst8(i):set8(i))+az).*rh528(sst8(i):set8(i))  /rholw*ms2mmd;
  
%   figure;
%   i=2;
%   plot(dt8d(sst8(i):set8(i))-dt8d(zp8at(i)),va(sst8(i):set8(i)),'b');hold on
%   plot(dt8d(sst8(i):set8(i))-dt8d(zp8at(i)),va1(sst8(i):set8(i)),'c');hold on

% water vapor flow due to advection
  advm(sst8(i):set8(i))=-va(sst8(i):set8(i)).*func.rhovs( (t528(sst8(i):set8(i))+t518(sst8(i):set8(i))+2*az)/2 ).*(rh528(sst8(i):set8(i)) +rh518(sst8(i):set8(i)) )/2/rholw;
  adv1(sst8(i):set8(i))=-va1(sst8(i):set8(i)).*func.rhovs( t518(sst8(i):set8(i))+az ).*rh518(sst8(i):set8(i)) /rholw;
  adv2(sst8(i):set8(i))=-va2(sst8(i):set8(i)).*func.rhovs( t528(sst8(i):set8(i))+az ).*rh528(sst8(i):set8(i)) /rholw;
  
  
  
%   pa=101325; %atmospheric pressure
%   tsk=40+az;
%   tak=60+az;
%   ptsk=func.svp(tsk);
%   ptak=func.svp(tak);
%   popf=pa*log( (pa-ptsk)/(pa-ptak))/(ptak-ptsk);
%   pf=(ptak-ptsk)/log( (pa-ptsk)/(pa-ptak));
%   
%        aaet=func.dv(  80+273.15    ).*...
%       (func.rhovs(80+az))...
%       /0.1/rholw*ms2mmd;

% the difference of water vapor between two sensors used for geting rac
% see page 203 notebook
drva(sst8(i):set8(i))=func.rhovs(t518(sst8(i):set8(i))+az).*rh518(sst8(i):set8(i))-func.rhovs(tma(1)+az).*rh(1);
% it is not necessary to use adv1 to obtain the changing aerodynamic 
% resistance
racc(sst8(i):set8(i))=drva(sst8(i):set8(i))/rholw./etsp8(sst8(i):set8(i));

uracc(sst8(i):set8(i))=(rh518(sst8(i):set8(i)).*func.drhovs(t518(sst8(i):set8(i))+az)*dtk+func.rhovs(t518(sst8(i):set8(i))+az)*drh )/rholw./etsp8(sst8(i):set8(i));

% uraccet assumes there are 10% difference due to evaporation. see page 207

uraccet(sst8(i):set8(i))=uracc(sst8(i):set8(i))-drva(sst8(i):set8(i))./rholw./etsp8(sst8(i):set8(i))*det;


% rsc(sst8(i):set8(i))=(func.svp(    t518(sst8(i):set8(i))+az     )*1./(  t518(sst8(i):set8(i))+az   )-func.svp(t528(sst8(i):set8(i)   )+az).*rh528(sst8(i):set8(i))./(t528(sst8(i):set8(i))+az)) ...
%       *mw/R/rholw./etsp8(sst8(i):set8(i))-racc(sst8(i):set8(i));  
rscc(sst8(i):set8(i))=(func.svp(    t518(sst8(i):set8(i))+az     )*1./(  t518(sst8(i):set8(i))+az   )-func.svp(tma(1)+az).*rh(1)./(tma(1)+az)) ...
      *mw/R/rholw./etsp8(sst8(i):set8(i))-racc(sst8(i):set8(i));    
  
  
  % there is almost no change to aerodynamic conditions comparing rsc and 
  % rsb
%     figure;
%     iv=10;
%         i=2; plot(mc8d(sst8(i):iv:set8(i)),racc(sst8(i):iv:set8(i)),'b+');hold on
%         plot(mc8d(sst8(i):iv:set8(i)),racc(sst8(i):iv:set8(i))+uracc(sst8(i):iv:set8(i)),'r*');hold on
        
        
        
        
%             figure;
%     iv=10;
%     i=4;    plot(mc8d(sst8(i):iv:set8(i)),uraccet(sst8(i):iv:set8(i))./racc(sst8(i):iv:set8(i)),'r*');hold on
        
        
        
        
%         i=3; plot(mc8d(sst8(i):iv:set8(i)),racc(sst8(i):iv:set8(i)),'ro');hold on
%         i=4; plot(mc8d(sst8(i):iv:set8(i)),racc(sst8(i):iv:set8(i)),'c*');hold on
%     
%     
%        
%        i=2; plot(mc8d(sst8(i):iv:set8(i)),racc(sst8(i):iv:set8(i)),'+');hold on
%        i=3; plot(mc8d(sst8(i):iv:set8(i)),racc(sst8(i):iv:set8(i)),'o');hold on
%        i=4; plot(mc8d(sst8(i):iv:set8(i)),racc(sst8(i):iv:set8(i)),'o');hold on
%        
%       figure;  
%       i=3;
%     plot(mc8d(sst8(i):set8(i)),rscc(sst8(i):set8(i)),'r'); hold on
%     plot(mc8d(sst8(i):set8(i)),rsb(sst8(i):set8(i)),'g'); hold on
%   
  
  
  
  end
  
  
  end
        fprintf(1, strcat('Interpolation finished \n'));

%% plot mean temperature against depth at each stable conditions

    figure
    semilogy(mc8d(sst8(2):iv:set8(2)),rsb(sst8(2):iv:set8(2),4),'go');hold on    %C9 140.16
    semilogy(mc8d(sst8(4):iv:set8(4)),rsb(sst8(4):iv:set8(4),4),'b.');   


%    h=figure;
%    set(gcf,'Units','normalized', 'WindowStyle','docked','OuterPosition',[0 0 1 1]); %maximized the window
%    subplot('Position',[0.05 0.65 0.18 0.32])
%    
% %  h1=plot([at52(1),at55(1),at51(1),at83(1),at84(1)],[es52,es55,es51,es83,es84],'r+','LineWidth',lw);hold on
%    i=2;h2=plot([at52(i),at55(i),at51(i),at83(i),at84(i)],[es52,es55,es51,es83,es84],'go','LineWidth',lw,'MarkerSize',ms);hold on
%    i=3;h3=plot([at52(i),at55(i),at51(i),at83(i),at84(i)],[es52,es55,es51,es83,es84],'b*','LineWidth',lw);hold on
% %   i=4;h4=plot([at52(i),at55(i),at51(i),at83(i),at84(i)],[es52,es55,es51,es83,es84],'cx','LineWidth',lw,'MarkerSize',ms);hold on
% %   i=5;h5=plot([at52(i),at55(i),at51(i),at83(i),at84(i)],[es52,es55,es51,es83,es84],'ms','LineWidth',lw,'MarkerSize',ms);hold on
% %   i=6;h6=plot([at52(i),at55(i),at51(i),at83(i),at84(i)],[es52,es55,es51,es83,es84],'yd','LineWidth',lw,'MarkerSize',ms);hold on
% %   i=7;h7=plot([at52(i),at55(i),at51(i),at83(i),at84(i)],[es52,es55,es51,es83,es84],'k^','LineWidth',lw,'MarkerSize',ms);hold on
% %   i=8;h8=plot([at52(i),at55(i),at51(i),at83(i),at84(i)],[es52,es55,es51,es83,es84],'rv','LineWidth',lw,'MarkerSize',ms);hold on
% %   i=9;h9=plot([at52(i),at55(i),at51(i),at83(i),at84(i)],[es52,es55,es51,es83,es84],'g<','LineWidth',lw,'MarkerSize',ms);hold on
% 
% 
% %   h8=plot([at5d(8),at5c(8),at5b(8),at5a(8)],[es5d,es5c,es5b,es5a],'rv','LineWidth',lw);hold on
% %   h9=plot([at5d(9),at5c(9),at5b(9),at5a(9)],[es5d,es5c,es5b,es5a],'g<','LineWidth',lw);hold on
% %   h10=plot([at5d(10),at5c(10),at5b(10),at5a(10)],[es5d,es5c,es5b,es5a],'b>','LineWidth',lw);hold on
%   
% %   f1=plot(ftp(1,1)*z1+ftp(1,2)-az,z1,'r:');hold on
% %   f2=plot(ftp(2,1)*z1+ftp(2,2)-az,z1,'g:');hold on
% %   f3=plot(ftp(3,1)*z1+ftp(3,2)-az,z1,'b:');hold on
% %   f4=plot(ftp(4,1)*z1+ftp(4,2)-az,z1,'c:');hold on
% %   f5=plot(ftp(5,1)*z1+ftp(5,2)-az,z1,'m:');hold on
% %   f6=plot(ftp(6,1)*z1+ftp(6,2)-az,z1,'y:');hold on
% %   f7=plot(ftp(7,1)*z1+ftp(7,2)-az,z1,'k:');hold on
% %   f8=plot(ftp(8,1)*z1+ftp(8,2)-az,z1,'r:');hold on
% %   f9=plot(ftp(9,1)*z1+ftp(9,2)-az,z1,'g:');hold on
% %   f10=plot(ftp(10,1)*z1+ftp(10,2)-az,z1,'b:');hold on
% %   
% %   f11=plot(fltp(1,1)*z2+fltp(1,2)-az,z2,'r:');hold on
% %   f12=plot(fltp(2,1)*z2+fltp(2,2)-az,z2,'g:');hold on
% %   f13=plot(fltp(3,1)*z2+fltp(3,2)-az,z2,'b:');hold on
% %   f14=plot(fltp(4,1)*z2+fltp(4,2)-az,z2,'c:');hold on
% %   f15=plot(fltp(5,1)*z2+fltp(5,2)-az,z2,'m:');hold on
% %   f16=plot(fltp(6,1)*z2+fltp(6,2)-az,z2,'y:');hold on
% %   f17=plot(fltp(7,1)*z2+fltp(7,2)-az,z2,'k:');hold on
% %   f18=plot(fltp(8,1)*z2+fltp(8,2)-az,z2,'r:');hold on
% %   f19=plot(fltp(9,1)*z2+fltp(9,2)-az,z2,'g:');hold on
% %   f20=plot(fltp(10,1)*z2+fltp(10,2)-az,z2,'b:');hold on
%   
%   
%   yl = get(gca,'YLim');
% %   line([atb(1) atb(1)],yl,'LineStyle',':','Color','-','LineWidth',lw);
% %   line([atb(2) atb(2)],yl,'LineStyle',':','Color','-','LineWidth',lw);
% %   line([atb(3) atb(3)],yl,'LineStyle',':','Color','b','LineWidth',lw);
% %   line([atb(4) atb(4)],yl,'LineStyle',':','Color','-','LineWidth',lw);
% %   line([atb(5) atb(5)],yl,'LineStyle',':','Color','-','LineWidth',lw);
% %   line([atb(6) atb(6)],yl,'LineStyle',':','Color','-','LineWidth',lw);
% %   line([atb(7) atb(7)],yl,'LineStyle',':','Color','-','LineWidth',lw);
% %   line([atb(8) atb(8)],yl,'LineStyle','-.','Color','-','LineWidth',lw);
% %   line([atb(9) atb(9)],yl,'LineStyle','-.','Color','-','LineWidth',lw);
% %   line([atb(10) atb(10)],yl,'LineStyle','-.','Color','b','LineWidth',lw);
%     
%  ax1=gca;set(ax1,'FontSize',tfs)  
%  xlabel('Temperature ( ^\circ C)','FontSize',fs,'FontWeight','bold')
%  ylabel('Depth (m)','FontSize',fs,'FontWeight','bold')
% %   hleg1 = legend([h2,h4,h5,h6],[cs(2),cs(4),cs(5),cs(6)],'Location','SouthEast');
% %   hleg1 = legend([h2,h3],[cs(2),cs(3)],'Location','SouthEast');
% %   set(hleg1, 'Box', 'on','FontSize',fs)
% % %% plot relative humidity of each sensors
% %    subplot('Position',[0.40 0.73 0.3 0.23])
% % %   h1=plot([arh5d(1),arh5c(1)],[es5d,es5c],'r+:','LineWidth',lw);hold on
% %   h2=plot([arh5d(2),arh5c(2)],[es5d,es5c],'go:','LineWidth',lw);hold on
% % %   h3=plot([arh5d(3),arh5c(3)],[es5d,es5c],'b*:','LineWidth',lw);hold on
% %   h4=plot([arh5d(4),arh5c(4)],[es5d,es5c],'cx:','LineWidth',lw);hold on
% %   h5=plot([arh5d(5),arh5c(5)],[es5d,es5c],'ms:','LineWidth',lw);hold on
% %   h6=plot([arh5d(6),arh5c(6)],[es5d,es5c],'yd:','LineWidth',lw);hold on
% % %   h7=plot([arh5d(7),arh5c(7)],[es5d,es5c],'k^:','LineWidth',lw);hold on
% % %   h8=plot([arh5d(8),arh5c(8)],[es5d,es5c],'rv:','LineWidth',lw);hold on
% % %   h9=plot([arh5d(9),arh5c(9)],[es5d,es5c],'g<:','LineWidth',lw);hold on
% % %   h10=plot([arh5d(10),arh5c(10)],[es5d,es5c],'b>:','LineWidth',lw);hold on
% %   yl = get(gca,'YLim');
% % %   line([arhb(1) arhb(1)],yl,'LineStyle',':','Color','-','LineWidth',lw);
% % %   line([arhb(2) arhb(2)],yl,'LineStyle',':','Color','-','LineWidth',lw);
% % %   line([arhb(3) arhb(3)],yl,'LineStyle',':','Color','b','LineWidth',lw);
% % %   line([arhb(4) arhb(4)],yl,'LineStyle',':','Color','-','LineWidth',lw);
% % %   line([arhb(5) arhb(5)],yl,'LineStyle',':','Color','-','LineWidth',lw);
% % %   line([arhb(6) arhb(6)],yl,'LineStyle',':','Color','-','LineWidth',lw);
% % %   line([arhb(7) arhb(7)],yl,'LineStyle',':','Color','-','LineWidth',lw);
% % %   line([arhb(8) arhb(8)],yl,'LineStyle','-.','Color','-','LineWidth',lw);
% % %   line([arhb(9) arhb(9)],yl,'LineStyle','-.','Color','-','LineWidth',lw);
% % %   line([arhb(10) arhb(10)],yl,'LineStyle','-.','Color','b','LineWidth',lw); 
% %   
% %  ax1=gca;set(ax1,'FontSize',tfs)  
% %  xlabel('Relative Humidity','FontSize',fs,'FontWeight','bold')
% %  ylabel('Depth (m)','FontSize',fs,'FontWeight','bold')
% % %  %% plot pressure of each sensors
% % %   subplot('Position',[0.40 0.06 0.3 0.45])
% % %   h1=plot([ap5b(1),ap5a(1),ap8b(1),ap8g(1)],[es5b,es5a,es8b,es8g],'r+:','LineWidth',lw);hold on
% % %   h2=plot([ap5b(2),ap5a(2),ap8b(2),ap8g(2)],[es5b,es5a,es8b,es8g],'go:','LineWidth',lw);hold on
% % %   h3=plot([ap5b(3),ap5a(3),ap8b(3),ap8g(3)],[es5b,es5a,es8b,es8g],'b*:','LineWidth',lw);hold on
% % %   h4=plot([ap5b(4),ap5a(4),ap8b(4),ap8g(4)],[es5b,es5a,es8b,es8g],'cx:','LineWidth',lw);hold on
% % %   h5=plot([ap5b(5),ap5a(5),ap8b(5),ap8g(5)],[es5b,es5a,es8b,es8g],'ms:','LineWidth',lw);hold on
% % %   h6=plot([ap5b(6),ap5a(6),ap8b(6),ap8g(6)],[es5b,es5a,es8b,es8g],'yd:','LineWidth',lw);hold on
% % %   h7=plot([ap5b(7),ap5a(7),ap8b(7),ap8g(7)],[es5b,es5a,es8b,es8g],'k^:','LineWidth',lw);hold on
% % %   h8=plot([ap5b(8),ap5a(8),ap8b(8),ap8g(8)],[es5b,es5a,es8b,es8g],'rv:','LineWidth',lw);hold on
% % %   h9=plot([ap5b(9),ap5a(9),ap8b(9),ap8g(9)],[es5b,es5a,es8b,es8g],'g<:','LineWidth',lw);hold on
% % %   h10=plot([ap5b(10),ap5a(10),ap8b(10),ap8g(10)],[es5b,es5a,es8b,es8g],'b>:','LineWidth',lw);hold on
% % %   
% % %     yl = get(gca,'YLim');
% % %  axis([-7 10 yl]) 
% % %     ax1=gca;set(ax1,'FontSize',tfs)  
% % %  xlabel('Pressure (m)','FontSize',fs,'FontWeight','bold')
% % %  ylabel('Depth (m)','FontSize',fs,'FontWeight','bold')
% 
% 
% %  %% plot moisture content of each sensors
% %   subplot('Position',[0.75 0.73 0.2 0.23])
% % %   h1=plot([amc8e(1),amc8f(1)],[es8e,es8f],'r+:','LineWidth',lw);hold on
% %   h2=plot([amc8e(2),amc8f(2)],[es8e,es8f],'go:','LineWidth',lw);hold on
% % %   h3=plot([amc8e(3),amc8f(3)],[es8e,es8f],'b*:','LineWidth',lw);hold on
% %   h4=plot([amc8e(4),amc8f(4)],[es8e,es8f],'cx:','LineWidth',lw);hold on
% %   h5=plot([amc8e(5),amc8f(5)],[es8e,es8f],'ms:','LineWidth',lw);hold on
% %   h6=plot([amc8e(6),amc8f(6)],[es8e,es8f],'yd:','LineWidth',lw);hold on
% % %   h7=plot([amc8e(7),amc8f(7)],[es8e,es8f],'k^:','LineWidth',lw);hold on
% % %   h8=plot([amc8e(8),amc8f(8)],[es8e,es8f],'rv:','LineWidth',lw);hold on
% % %   h9=plot([amc8e(9),amc8f(9)],[es8e,es8f],'g<:','LineWidth',lw);hold on
% %   h10=plot([amc8e(10),amc8f(10)],[es8e,es8f],'b>:','LineWidth',lw);hold on
% %   yl = get(gca,'YLim');
% %    ax1=gca;set(ax1,'FontSize',tfs)  
% %  xlabel('Saturation','FontSize',fs,'FontWeight','bold')
% %  ylabel('Depth (m)','FontSize',fs,'FontWeight','bold')
%  saveas(h,['Coarse_90_Analysis_Transient_Temperature','fig'],'fig')
%   fprintf(1, strcat('First Graph finished \n'));
% %% plot Moisture content after translation
%  h=figure;
%  set(gcf,'Units','normalized', 'WindowStyle','docked','OuterPosition',[0 0 1 1]); %maximized the window
% 
% 
% 
% subplot('Position',[0.55 0.70 0.43 0.28])
% for i=1:ns
%  if (wl(i)==0 && (sc(i)==2 ||sc(i)==3))
%   plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),mc8d(sst8(i):set8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
%  end
% end
% 
% %i=2;h2=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),mc8d(zp8at(i):set8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% %i=3;h3=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),mc8d(zp8at(i):set8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% %i=4;h4=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),mc8d(zp8at(i):set8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=5;h5=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),mc8d(zp8at(i):set8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=6;h6=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),mc8d(zp8at(i):set8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=7;h7=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),mc8d(zp8at(i):set8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=8;h8=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),mc8d(zp8at(i):set8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=8;h9=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),mc8d(zp8at(i):set8(i)),'g--','LineWidth',lw,'color',col(i,:));hold on
% 
% axis([0,8,0,1])
% xlabel('Time (day)','FontSize',fs)
% ylabel('Sautration (-)','FontSize',fs)
% %% plot accumulative evaporation rate after translation
% subplot('Position',[0.55 0.38 0.43 0.28])
% for i=1:ns
%  if (wl(i)==0 && (sc(i)==2 ||sc(i)==3))
% %plot(dtsd(ssts(i):sets(i))-dtsd(ssts(i)),aet9c(ssts(i):sets(i))-aet9c(ssts(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),aetsp8(sst8(i):set8(i))-aetsp8(sst8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
%  end
% end
% %i=2;h2=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% %i=3;h3=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=4;h4=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=5;h5=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=6;h6=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=7;h7=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=8;h8=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=9;h9=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),'g--','LineWidth',lw,'color',col(i,:));hold on
% 
% ylabel('AET (mm)','FontSize',fs)
% xlabel('Time (day)','FontSize',fs)
% %% Plot Transient evaportion rate after translation as a function of time
% subplot('Position',[0.55 0.05 0.43 0.28])
% 
% for i=1:ns
%  if (wl(i)==0 && (sc(i)==2 ||sc(i)==3))
% % plot(dtsd(ssts(i):sets(i))-dtsd(ssts(i)),etsgf9(ssts(i):sets(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),etsp8(sst8(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% %plot(dt8d(sst8(i):iv:set8(i))-dt8d(sst8(i)),adv1(sst8(i):iv:set8(i)),'*','LineWidth',lw,'color',col(i,:));hold on
% %plot(dt8d(sst8(i):iv:set8(i))-dt8d(sst8(i)),adv2(sst8(i):iv:set8(i)),'o','LineWidth',lw,'color',col(i,:));hold on
% plot(dt8d(sst8(i):iv:set8(i))-dt8d(sst8(i)),adv1(sst8(i):iv:set8(i))*ms2mmd,'+','LineWidth',lw,'color',col(i,:));hold on
%  end
% end
% %i=2;h2=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf9(zpsat(i):sets(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% %i=3;h3=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf9(zpsat(i):sets(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=4;h4=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf9(zpsat(i):sets(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=5;h5=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf9(zpsat(i):sets(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=6;h6=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf9(zpsat(i):sets(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=7;h7=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf9(zpsat(i):sets(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=8;h8=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf9(zpsat(i):sets(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=9;h9=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf9(zpsat(i):sets(i))*ms2mmd,'g--','LineWidth',lw,'color',col(i,:));hold on
% 
% xlabel('Time (day)','FontSize',fs)
% ylabel('ET rate (mm/day)','FontSize',fs)
% %% Plot Transient evaportion rate as a function of accumulative evaporation rate
% subplot('Position',[0.04 0.68 0.43 0.28])
% 
% for i=1:ns
%  if (wl(i)==0 && (sc(i)==2 ||sc(i)==3))
% plot(aetsp8(sst8(i):set8(i)),etsp8(sst8(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
%  end
% end
% 
% %i=2;h2=plot(aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),etsgf9(zpsat(i):sets(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% %i=3;h3=plot(aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),etsgf9(zpsat(i):sets(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=4;h4=plot(aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),etsgf9(zpsat(i):sets(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=5;h5=plot(aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),etsgf9(zpsat(i):sets(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=6;h6=plot(aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),etsgf9(zpsat(i):sets(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=7;h7=plot(aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),etsgf9(zpsat(i):sets(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=8;h8=plot(aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),etsgf9(zpsat(i):sets(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=9;h9=plot(aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),etsgf9(zpsat(i):sets(i))*ms2mmd,'g--','LineWidth',lw,'color',col(i,:));hold on
% 
% xlabel('AET (mm)','FontSize',fs)
% ylabel('ET rate (mm/day)','FontSize',fs)
% %% Plot Transient evaporation rate as a function of saturation
% subplot('Position',[0.04 0.36 0.43 0.28])
% 
% for i=1:ns
%  if (wl(i)==0 && (sc(i)==2 ||sc(i)==3))
% plot(mc8d(sst8(i):set8(i)),etsp8(sst8(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
%  end
% end
% 
% 
% %i=2;h2=plot(mc8d(zp8at(i):set8(i)),etsp8(zp8at(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% %i=3;h3=plot(mc8d(zp8at(i):set8(i)),etsp8(zp8at(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=4;h4=plot(mc8d(zp8at(i):set8(i)),etsp8(zp8at(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=5;h5=plot(mc8d(zp8at(i):set8(i)),etsp8(zp8at(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=6;h6=plot(mc8d(zp8at(i):set8(i)),etsp8(zp8at(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=7;h7=plot(mc8d(zp8at(i):set8(i)),etsp8(zp8at(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=8;h8=plot(mc8d(zp8at(i):set8(i)),etsp8(zp8at(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=9;h9=plot(mc8d(zp8at(i):set8(i)),etsp8(zp8at(i):set8(i))*ms2mmd,'g--','LineWidth',lw,'color',col(i,:));hold on
% 
% xlabel('Saturation (-)','FontSize',fs)
% ylabel('ET rate (mm/day)','FontSize',fs)
% 
% 
% 
% 
% %% Plot Transient evaporation rate based on the relative humidity
% subplot('Position',[0.04 0.04 0.43 0.28])
% 
% for i=1:ns
%  if (wl(i)==0 && (sc(i)==2 ||sc(i)==3))
% plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),etrh8(sst8(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% 
% %plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),etrh812(sst8(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% 
% 
%  end
% end
% 
% %i=2;h2=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),etrh8(zp8at(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% %i=3;h3=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),etrh8(zp8at(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=4;h4=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),etrh8(zp8at(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=5;h5=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),etrh8(zp8at(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=6;h6=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),etrh8(zp8at(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=7;h7=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),etrh8(zp8at(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=8;h8=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),etrh8(zp8at(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% % i=9;h9=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),etrh8(zp8at(i):set8(i))*ms2mmd,'g--','LineWidth',lw,'color',col(i,:));hold on
% 
% xlabel('Time (day)','FontSize',fs)
% ylabel('ET by RH (mm/day)','FontSize',fs)
% axis([0,8,0,1])
% 
% 
% saveas(h,['Coarse_90_Analysis_Transient_ETrate','.fig'],'fig')
% fprintf(1, strcat('Second Graph finished \n'));
% 
% 
% 
% % Plot the non-dimensional result
% %% Plot non-dimensional evaporation rate as a function of accumulative evaporation 
% h=figure;
% set(gcf,'Units','normalized', 'WindowStyle','docked','OuterPosition',[0 0 1 1]); %maximized the window
% subplot('Position',[0.04 0.70 0.43 0.28])
% 
% for i=1:ns
%  if (wl(i)==0 && (sc(i)==2 ||sc(i)==3))
% plot(aetsp8(sst8(i):set8(i)),etsp8(sst8(i):set8(i))/etsp8(sst8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
%  end
% end
% 
% 
% 
% %i=2;h2=plot(aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),etsgf9(zpsat(i):sets(i))/etsgf9(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% %i=3;h3=plot(aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),etsgf9(zpsat(i):sets(i))/etsgf9(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=4;h4=plot(aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),etsgf9(zpsat(i):sets(i))/etsgf9(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=5;h5=plot(aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),etsgf9(zpsat(i):sets(i))/etsgf9(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=6;h6=plot(aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),etsgf9(zpsat(i):sets(i))/etsgf9(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=7;h7=plot(aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),etsgf9(zpsat(i):sets(i))/etsgf9(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=8;h8=plot(aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),etsgf9(zpsat(i):sets(i))/etsgf9(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=9;h9=plot(aet9c(zpsat(i):sets(i))-aet9c(zpsat(i)),etsgf9(zpsat(i):sets(i))/etsgf9(zpsat(i)),'g--','LineWidth',lw,'color',col(i,:));hold on
% 
% xlabel('AET (mm)','FontSize',fs)
% ylabel('ET rate (-)','FontSize',fs)
% 
% 
% %% Plot non-dimensional evaporation rate as a function of saturation
% 
% subplot('Position',[0.04 0.36 0.43 0.28])
% for i=1:ns
%  if (wl(i)==0 && (sc(i)==2 ||sc(i)==3))
% plot(mc8d(sst8(i):set8(i)),etsp8(sst8(i):set8(i))/etsp8(sst8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
%  end
% end
% xlabel('Saturation (-)','FontSize',fs)
% ylabel('ET rate (-)','FontSize',fs)
% axis([0,1,0,1])
% %% Plot non-dimensional evaporation rate calculated by relative humidity against time
% subplot('Position',[0.04 0.04 0.43 0.28])
% 
% 
% j=0;
% for i=1:ns
%  if (wl(i)==0 && (sc(i)==2 ||sc(i)==3))
% 
% plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),etrh8(sst8(i):set8(i))/etrh8(sst8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% 
% 
% plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),etrh812(sst8(i):set8(i))/etrh812(sst8(i)),'b','LineWidth',lw,'color',col(i,:));hold on
%  end
% end
% 
% %i=2;h2=plot(dt8d(zp8at(i)+j:set8(i))-dt8d(zp8at(i)+j),etrh8(zp8at(i)+j:set8(i))/etrh8(zp8at(i)+j),'-','LineWidth',lw,'color',col(i,:));hold on
% %i=3;h3=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),etrh8(zp8at(i):set8(i))/etrh8(zp8at(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=4;h4=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),etrh8(zp8at(i):set8(i))/etrh8(zp8at(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=5;h5=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),etrh8(zp8at(i):set8(i))/etrh8(zp8at(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=6;h6=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),etrh8(zp8at(i):set8(i))/etrh8(zp8at(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=7;h7=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),etrh8(zp8at(i):set8(i))/etrh8(zp8at(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=8;h8=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),etrh8(zp8at(i):set8(i))/etrh8(zp8at(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=9;h9=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),etrh8(zp8at(i):set8(i))/etrh8(zp8at(i)),'g--','LineWidth',lw,'color',col(i,:));hold on
% 
% xlabel('Time (day)','FontSize',fs)
% ylabel('ET by RH (-)','FontSize',fs)
% axis([0,8,0,1])
% 
% 
% %% Plot non-dimensional evaporation rate calculated by relative humidity against time
% subplot('Position',[0.55 0.36 0.43 0.28])
% 
% 
% for i=1:ns
%  if (wl(i)==0 && (sc(i)==2 ||sc(i)==3))
% plot(aetsp8(sst8(i):set8(i)),etrh8(sst8(i):set8(i))/etrh8(sst8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% 
%  end
% end
% 
% %i=2;h2=plot(aet9c8(zp8at(i):set8(i))-aet9c8(zp8at(i)),etrh8(zp8at(i):set8(i))/etrh8(zp8at(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% %i=3;h3=plot(aet9c8(zp8at(i):set8(i))-aet9c8(zp8at(i)),etrh8(zp8at(i):set8(i))/etrh8(zp8at(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=4;h4=plot(aet9c8(zp8at(i):set8(i))-aet9c8(zp8at(i)),etrh8(zp8at(i):set8(i))/etrh8(zp8at(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=5;h5=plot(aet9c8(zp8at(i):set8(i))-aet9c8(zp8at(i)),etrh8(zp8at(i):set8(i))/etrh8(zp8at(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=6;h6=plot(aet9c8(zp8at(i):set8(i))-aet9c8(zp8at(i)),etrh8(zp8at(i):set8(i))/etrh8(zp8at(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=7;h7=plot(aet9c8(zp8at(i):set8(i))-aet9c8(zp8at(i)),etrh8(zp8at(i):set8(i))/etrh8(zp8at(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=8;h8=plot(aet9c8(zp8at(i):set8(i))-aet9c8(zp8at(i)),etrh8(zp8at(i):set8(i))/etrh8(zp8at(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=9;h9=plot(aet9c8(zp8at(i):set8(i))-aet9c8(zp8at(i)),etrh8(zp8at(i):set8(i))/etrh8(zp8at(i)),'g--','LineWidth',lw,'color',col(i,:));hold on
% xlabel('AET (mm)','FontSize',fs)
% ylabel('ET by RH (-)','FontSize',fs)
% 
% 
% %% Plot Transient evaportion rate after translation as a function of time
% subplot('Position',[0.55 0.05 0.43 0.28])
% 
% for i=1:ns
%  if (wl(i)==0 && (sc(i)==2 ||sc(i)==3))
% plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),etsp8(sst8(i):set8(i))/etsp8(sst8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% plot(dt8d(sst8(i):iv:set8(i))-dt8d(sst8(i)),adv1(sst8(i):iv:set8(i))/adv1(sst8(i)),'*','LineWidth',lw,'color',col(i,:));hold on
% plot(dt8d(sst8(i):iv:set8(i))-dt8d(sst8(i)),adv2(sst8(i):iv:set8(i))/adv2(sst8(i)),'o','LineWidth',lw,'color',col(i,:));hold on
% plot(dt8d(sst8(i):iv:set8(i))-dt8d(sst8(i)),advm(sst8(i):iv:set8(i))/advm(sst8(i)),'+','LineWidth',lw,'color',col(i,:));hold on
%  end
% end
% 
% 
% %i=2;h2=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf9(zpsat(i):sets(i))/etsgf9(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% %i=3;h3=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf9(zpsat(i):sets(i))/etsgf9(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=4;h4=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf9(zpsat(i):sets(i))/etsgf9(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=5;h5=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf9(zpsat(i):sets(i))/etsgf9(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=6;h6=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf9(zpsat(i):sets(i))/etsgf9(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=7;h7=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf9(zpsat(i):sets(i))/etsgf9(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=8;h8=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf9(zpsat(i):sets(i))/etsgf9(zpsat(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % i=9;h9=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf9(zpsat(i):sets(i))/etsgf9(zpsat(i)),'g--','LineWidth',lw,'color',col(i,:));hold on
% 
% xlabel('Time (day)','FontSize',fs)
% ylabel('non-ET rate (-)','FontSize',fs)
% 
% 
% saveas(h,['Coarse_90_Analysis_Transient_nondimensional','.fig'],'fig')
% fprintf(1, strcat('non-dimensional Graph finished \n'));
% % %% Plot Surface resistance as a function of saturation
% % h=figure;
% % set(gcf,'Units','normalized', 'WindowStyle','docked','OuterPosition',[0 0 1 1]); %maximized the window
% % 
% % %  subplot('Position',[0.04 0.70 0.43 0.28])
% % %   ha=plot(0.01:0.05:1,func.rs1994(0.01:0.05:1,0.375,por),'k-','LineWidth',lw,'color',col(i,:));hold on % van der grivend 1994
% % %   hb=plot(0.03:0.01:0.45,func.rs1986(0.03:0.01:0.45,por),'k-.','LineWidth',lw,'color',col(i,:));hold on % camilo 1986
% % %   hc=plot(0.01:0.05:1,func.rs1996(0.01:0.05:1,por),':','LineWidth',lw,'color',col(i,:));hold on % Daamen 1996
% % %   hd=plot(0.01:0.05:1,func.rs1984(0.01:0.05:1),'--','LineWidth',lw,'color',col(i,:));hold on % Sun 1984
% % %   
% % %   
% % %   i=2;h2=plot(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % %   i=4;h4=plot(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % %   i=5;h5=plot(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % %   i=6;h6=plot(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % %   i=7;h7=plot(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % %   i=8;h8=plot(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% % %   
% % %   
% % % 
% % %   
% % %   axis([0,1,0,600])
% % %   hleg1 = legend([ha,hb,hc,hd],{sprintf('Van der Griend et al. [1994]'),sprintf('Camillo et al. [1986]'),sprintf('Daamen et al. [1996]'), ...
% % %     sprintf('Sun [1984]')},'Location','NorthEast');
% % %   xlabel('Saturation (-)','FontSize',fs)
% % %   ylabel('rs (s/m)','FontSize',fs)
% % 
% % %% Plot surface resistance on log profile
% % subplot('Position',[0.04 0.55 0.45 0.4])
% %   ha=semilogy(0.01:0.05:1,func.rs1994(0.01:0.05:1,0.375,por),'k-','LineWidth',lw,'color',col(i,:));hold on % van der grivend 1994
% %   hb=semilogy(0.03:0.01:0.45,func.rs1986(0.03:0.01:0.45,por),'k-.','LineWidth',lw,'color',col(i,:));hold on % camilo 1986
% %   hc=semilogy(0.01:0.05:1,func.rs1996(0.01:0.05:1,por),':','LineWidth',lw,'color',col(i,:));hold on % Daamen 1996
% %   hd=semilogy(0.01:0.05:1,func.rs1984(0.01:0.05:1),'--','LineWidth',lw,'color',col(i,:));hold on % Sun 1984
% %   
% %   
% %   j=2;
% %   i=2;h2=semilogy(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i),j),'-','LineWidth',lw,'color',col(i,:)2);hold on
% %   i=3;h3=semilogy(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i),j),'-','LineWidth',lw,'color',col(i,:)2);hold on  
% % %   i=4;h4=semilogy(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i),j),'-','LineWidth',lw,'color',col(i,:)2);hold on
% % %   i=5;h5=semilogy(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i),j),'-','LineWidth',lw,'color',col(i,:)2);hold on
% % %   i=6;h6=semilogy(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i),j),'-','LineWidth',lw,'color',col(i,:)2);hold on
% % %   i=7;h7=semilogy(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i),j),'-','LineWidth',lw,'color',col(i,:)2);hold on
% % %   i=8;h8=semilogy(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i),j),'-','LineWidth',lw,'color',col(i,:)2);hold on
% % %   i=9;h9=semilogy(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i),j),'g--','LineWidth',lw,'color',col(i,:)2);hold on
% %   
% %   i=2;g22=plot(mc8d(hsp8(i):set8(i)),rshs2(hsp8(i):set8(i)),'g--','LineWidth',lw,'color',col(i,:)2);hold on
% %   i=3;g32=plot(mc8d(hsp8(i):set8(i)),rshs2(hsp8(i):set8(i)),'c--','LineWidth',lw,'color',col(i,:)2);hold on  
% % %   i=4;g42=plot(mc8d(hsp8(i):set8(i)),rshs2(hsp8(i):set8(i)),'c--','LineWidth',lw,'color',col(i,:)2);hold on
% % %   i=5;g52=plot(mc8d(hsp8(i):set8(i)),rshs2(hsp8(i):set8(i)),'m--','LineWidth',lw,'color',col(i,:)2);hold on
% % %   i=6;g62=plot(mc8d(hsp8(i):set8(i)),rshs2(hsp8(i):set8(i)),'y--','LineWidth',lw,'color',col(i,:)2);hold on
% % %   i=7;g72=plot(mc8d(hsp8(i):set8(i)),rshs2(hsp8(i):set8(i)),'k--','LineWidth',lw,'color',col(i,:)2);hold on
% % %   i=8;g82=plot(mc8d(hsp8(i):set8(i)),rshs2(hsp8(i):set8(i)),'r--','LineWidth',lw,'color',col(i,:)2);hold on
% % %   i=9;g92=plot(mc8d(hsp8(i):set8(i)),rshs2(hsp8(i):set8(i)),'g--','LineWidth',lw,'color',col(i,:)2);hold on
% %   
% %   
% %   axis([0,1,1,5000])
% %   hleg1 = legend([ha,hb,hc,hd],{sprintf('Van der Griend et al. [1994]'),sprintf('Camillo et al. [1986]'),sprintf('Daamen et al. [1996]'), ...
% %   sprintf('Sun [1984]')},'Location','NorthEast');
% %   xlabel('Saturation (-)','FontSize',fs)
% %   ylabel('rs (s/m)','FontSize',fs)
% % %% Plot Surface resistance as a function of saturation with translation at the full saturation
% %  subplot('Position',[0.04 0.04 0.45 0.4])
% %   ha=plot(0.01:0.05:1,func.rs1994(0.01:0.05:1,0.375,por),'k-','LineWidth',lw,'color',col(i,:));hold on % van der grivend 1994
% %   hb=plot(0.03:0.01:0.45,func.rs1986(0.03:0.01:0.45,por),'k-.','LineWidth',lw,'color',col(i,:));hold on % camilo 1986
% %   hc=plot(0.01:0.05:1,func.rs1996(0.01:0.05:1,por),':','LineWidth',lw,'color',col(i,:));hold on % Daamen 1996
% %   hd=plot(0.01:0.05:1,func.rs1984(0.01:0.05:1),'--','LineWidth',lw,'color',col(i,:));hold on % Sun 1984
% %   
% % %  surface resistance from half saturation  
% %   i=2;g2=plot(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'g-.','LineWidth',lw,'color',col(i,:)2);hold on
% %   i=3;g3=plot(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'c-.','LineWidth',lw,'color',col(i,:)2);hold on
% % %   i=4;g4=plot(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'c-.','LineWidth',lw,'color',col(i,:)2);hold on  
% % %   i=5;g5=plot(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'m-.','LineWidth',lw,'color',col(i,:)2);hold on
% % %   i=6;g6=plot(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'y-.','LineWidth',lw,'color',col(i,:)2);hold on
% % %   i=7;g7=plot(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'k-.','LineWidth',lw,'color',col(i,:)2);hold on
% % %   i=8;g8=plot(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'r-.','LineWidth',lw,'color',col(i,:)2);hold on
% % %   i=9;g9=plot(mc8d(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'g-.','LineWidth',lw,'color',col(i,:)2);hold on
% %   
% %   i=2;h2=plot(mc8d(zp8at(i):set8(i)),rs2(zp8at(i):set8(i)),'-','LineWidth',lw2);hold on
% %   i=3;h3=plot(mc8d(zp8at(i):set8(i)),rs2(zp8at(i):set8(i)),'-','LineWidth',lw2);hold on  
% % %   i=4;h4=plot(mc8d(zp8at(i):set8(i)),rs2(zp8at(i):set8(i)),'-','LineWidth',lw2);hold on
% % %   i=5;h5=plot(mc8d(zp8at(i):set8(i)),rs2(zp8at(i):set8(i)),'-','LineWidth',lw2);hold on
% % %   i=6;h6=plot(mc8d(zp8at(i):set8(i)),rs2(zp8at(i):set8(i)),'-','LineWidth',lw2);hold on
% % %   i=7;h7=plot(mc8d(zp8at(i):set8(i)),rs2(zp8at(i):set8(i)),'-','LineWidth',lw2);hold on
% % %   i=8;h8=plot(mc8d(zp8at(i):set8(i)),rs2(zp8at(i):set8(i)),'-','LineWidth',lw2);hold on
% % %   i=9;h9=plot(mc8d(zp8at(i):set8(i)),rs2(zp8at(i):set8(i)),'-','LineWidth',lw2);hold on
% % 
% %    
% %   i=2;hh2=plot(mc8d(zp8at(i):set8(i)),rs3(zp8at(i):set8(i)),'g--','LineWidth',lw2);hold on
% %   i=3;hh3=plot(mc8d(zp8at(i):set8(i)),rs3(zp8at(i):set8(i)),'c--','LineWidth',lw2);hold on
% % %   i=4;hh4=plot(mc8d(zp8at(i):set8(i)),rs3(zp8at(i):set8(i)),'c--','LineWidth',lw2);hold on
% % %   i=5;hh5=plot(mc8d(zp8at(i):set8(i)),rs3(zp8at(i):set8(i)),'m--','LineWidth',lw2);hold on
% % %   i=6;hh6=plot(mc8d(zp8at(i):set8(i)),rs3(zp8at(i):set8(i)),'y--','LineWidth',lw2);hold on
% % %   i=7;hh7=plot(mc8d(zp8at(i):set8(i)),rs3(zp8at(i):set8(i)),'k--','LineWidth',lw2);hold on
% % %   i=8;hh8=plot(mc8d(zp8at(i):set8(i)),rs3(zp8at(i):set8(i)),'r--','LineWidth',lw2);hold on
% % %   i=9;hh9=plot(mc8d(zp8at(i):set8(i)),rs3(zp8at(i):set8(i)),'g--','LineWidth',lw2);hold on
% %   
% %   
% %   axis([0,1,-100,600])
% %   hleg1 = legend([ha,hb,hc,hd,g2,h2,hh2],{sprintf('Van der Griend et al. [1994]'),sprintf('Camillo et al. [1986]'),sprintf('Daamen et al. [1996]'), ...
% %   sprintf('Sun [1984]'),sprintf('Tma=%3.2f, hra=%3.2f',tma(1),rh(1)),sprintf('Tma=%3.2f, hra=%3.2f',tma(2),rh(2)), ...
% %   sprintf('Tma=%3.2f, hra=%3.2f',tma(3),rh(3))},'Location','NorthEast');
% %   xlabel('Saturation (-)','FontSize',fs)
% %   ylabel('rs (s/m)','FontSize',fs)
% %   
% %   
%   %% Plot Surface resistance as a function of saturation with translation at the half saturation
% %   figure;
% %   %  surface resistance from half saturation
% %   % subplot('Position',[0.54 0.04 0.45 0.4])
% %   ha=plot(0.01:0.05:1,func.rs1994(0.01:0.05:1,0.375,por),'k-','LineWidth',lw);hold on % van der grivend 1994
% %   hb=plot(0.03:0.01:0.45,func.rs1986(0.03:0.01:0.45,por),'k-.','LineWidth',lw);hold on % camilo 1986
% %   hc=plot(0.01:0.05:1,func.rs1996(0.01:0.05:1,por),':','LineWidth',lw);hold on % Daamen 1996
% %   hd=plot(0.01:0.05:1,func.rs1984(0.01:0.05:1),'--','LineWidth',lw);hold on % Sun 1984
% %   
% %   i=4;h3=plot(mc8d(sst8(i):iv:set8(i)),rsb(sst8(i):iv:set8(i),l),'.','LineWidth',lw,'color',col(i,:));hold on
% % k=1;h3=plot(mc8d(zp8at(i):iv:set8(i)),rscc(zp8at(i):iv:set8(i),k)-rscc(zp8at(i),k),'--','LineWidth',lw,'color',col(i,:));hold on
% % 
% % i=2;h3=plot(mc8d(sst8(i):iv:set8(i)),rsb(sst8(i):iv:set8(i),l),'.','LineWidth',lw,'color',col(i,:));hold on
% % k=2;h3=plot(mc8d(zp8at(i):iv:set8(i)),rscc(zp8at(i):iv:set8(i),k)-rscc(zp8at(i),k),'-.','LineWidth',lw,'color',col(i,:));hold on
% % 
% %   axis([0,1,0,1000])
% %    
% %    
% %   i=2;g2=plot(mc8d(hsp8(i):set8(i)),rshs(hsp8(i):set8(i)),'g-.','LineWidth',lw2);hold on
% %   i=3;g3=plot(mc8d(hsp8(i):set8(i)),rshs(hsp8(i):set8(i)),'c-.','LineWidth',lw2);hold on
% % %   i=4;g4=plot(mc8d(hsp8(i):set8(i)),rshs(hsp8(i):set8(i)),'c-.','LineWidth',lw2);hold on
% % %   i=5;g5=plot(mc8d(hsp8(i):set8(i)),rshs(hsp8(i):set8(i)),'m-.','LineWidth',lw2);hold on
% % %   i=6;g6=plot(mc8d(hsp8(i):set8(i)),rshs(hsp8(i):set8(i)),'y-.','LineWidth',lw2);hold on
% % %   i=7;g7=plot(mc8d(hsp8(i):set8(i)),rshs(hsp8(i):set8(i)),'k-.','LineWidth',lw2);hold on
% % %   i=8;g8=plot(mc8d(hsp8(i):set8(i)),rshs(hsp8(i):set8(i)),'r-.','LineWidth',lw2);hold on
% % %   i=9;g9=plot(mc8d(hsp8(i):set8(i)),rshs(hsp8(i):set8(i)),'g-.','LineWidth',lw2);hold on
% %   
% %   
% %   i=2;h2=plot(mc8d(hsp8(i):set8(i)),rshs2(hsp8(i):set8(i)),'-','LineWidth',lw2);hold on
% %   i=3;h3=plot(mc8d(hsp8(i):set8(i)),rshs2(hsp8(i):set8(i)),'-','LineWidth',lw2);hold on  
% % %   i=4;h4=plot(mc8d(hsp8(i):set8(i)),rshs2(hsp8(i):set8(i)),'-','LineWidth',lw2);hold on
% % %   i=5;h5=plot(mc8d(hsp8(i):set8(i)),rshs2(hsp8(i):set8(i)),'-','LineWidth',lw2);hold on
% % %   i=6;h6=plot(mc8d(hsp8(i):set8(i)),rshs2(hsp8(i):set8(i)),'-','LineWidth',lw2);hold on
% % %   i=7;h7=plot(mc8d(hsp8(i):set8(i)),rshs2(hsp8(i):set8(i)),'-','LineWidth',lw2);hold on
% % %   i=8;h8=plot(mc8d(hsp8(i):set8(i)),rshs2(hsp8(i):set8(i)),'-','LineWidth',lw2);hold on
% % %   i=9;h9=plot(mc8d(hsp8(i):set8(i)),rshs2(hsp8(i):set8(i)),'-','LineWidth',lw2);hold on
% % %   
% %   i=2;hh2=plot(mc8d(hsp8(i):set8(i)),rshs3(hsp8(i):set8(i)),'g--','LineWidth',lw2);hold on
% %   i=3;hh3=plot(mc8d(hsp8(i):set8(i)),rshs3(hsp8(i):set8(i)),'c--','LineWidth',lw2);hold on  
% % %   i=4;hh4=plot(mc8d(hsp8(i):set8(i)),rshs3(hsp8(i):set8(i)),'c--','LineWidth',lw2);hold on
% % %   i=5;hh5=plot(mc8d(hsp8(i):set8(i)),rshs3(hsp8(i):set8(i)),'m--','LineWidth',lw2);hold on
% % %   i=6;hh6=plot(mc8d(hsp8(i):set8(i)),rshs3(hsp8(i):set8(i)),'y--','LineWidth',lw2);hold on
% % %   i=7;hh7=plot(mc8d(hsp8(i):set8(i)),rshs3(hsp8(i):set8(i)),'k--','LineWidth',lw2);hold on
% % %   i=8;hh8=plot(mc8d(hsp8(i):set8(i)),rshs3(hsp8(i):set8(i)),'r--','LineWidth',lw2);hold on
% % %   i=9;hh9=plot(mc8d(hsp8(i):set8(i)),rshs3(hsp8(i):set8(i)),'g--','LineWidth',lw2);hold on
% %   axis([0,1,-100,600])
% %   hleg1 = legend([ha,hb,hc,hd,g2,h2,hh2],{sprintf('Van der Griend et al. [1994]'),sprintf('Camillo et al. [1986]'),sprintf('Daamen et al. [1996]'), ...
% %   sprintf('Sun [1984]'),sprintf('Tma=%3.2f, hra=%3.2f',tma(1),rh(1)),sprintf('Tma=%3.2f, hra=%3.2f',tma(2),rh(2)), ...
% %   sprintf('Tma=%3.2f, hra=%3.2f',tma(3),rh(3))},'Location','NorthEast');
% %   xlabel('Saturation (-)','FontSize',fs)
% %   ylabel('rs (s/m)','FontSize',fs)
% %   
% %  %% save file
% %  saveas(h,['Coarse_90_Analysis_Transient_Resistance','.fig'],'fig')
% save('Coarse.mat','-v4') 
%  fprintf(1, strcat('Analysis finished \n'));

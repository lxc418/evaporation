%  Comparision against the different senarios
lw=2;tfs=10;fs=12;ms=5;por=0.39;
tma=[60;80;100];rh=[0.05;0.01;0];mw=0.018;R=8.314;
az=273.15;
kg=5e-12; % permeability of air (m2) by thomas
miua=1.846e-5; % viscosity of air (kg/m/s) by thomas
dz=0.3;    % meter
dz1=0.2;
iv=10;
%% get the elevation of each sensor
% (E)levation of (S)ensors from DT(8)5g at channel (C) (unit m) for moisture content
% (E)levation of (S)ensors from DT(8)5g at channel (1) (unit m) for temperature and pore water pressure
   es8c=-0.01;es8d=-0.05;es8e=-0.01;es8f=-0.05;es81=-0.08;es82=-0.03;es83=-0.03;es84=-0.08;
   % (E)levation of (S)ensors from EM(5)0 at Channel (A) 1-a 2-b 3-c 4-d ... (unit m)
   es51=0.01;es52=0.05;es53=0.01;es54=0.02;es55=0.02;
   % (E)levation at the (B)ottom of the (9)0mm (C)olumn
   eb9c=-0.2;
   % elevation of the light
% ALL THE UNIT SHOULD BE IN STANDARD UNIT (M KG K)
% Change all of unit of the average temperature values from celsius to kelvin
% at5ak=at54+az;at5bk=at53+az;
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
%% calculate surface resistance
% the time series of surface resistance is based on the dt85g time series
% (A)ccumulative (E)vapora(T)ion of (6)0 mm (C)olumn rate with respected to scale data converted to DT(8)5G time data 
 aet6c8=zeros(nr8,1);
% (T)ransient (E)vapora(T)ion of 60 mm column rate with respected to (S)cale data converted to(2) DT(8)5G time data 
 et68=zeros(nr8,1);
% (E)vapora(T)ion rate at (6)0mm column calculated through the (R)elative (H)umidity sensor at the tunnel using time series of DT(8)5G
 et6rh8=zeros(nr8,1);
% (R)elative (H)umidity at EM(5)0 Channel (3) interpolated into the time series of DT(8)5G 
 rh538=zeros(nr8,1);
 rh548=zeros(nr8,1);
% (T)emperature at EM(5)0 Channel (3) interpolated into the time series of DT(8)5G 
 t538=zeros(nr8,1);
 t548=zeros(nr8,1);
 
% (RS) using (T)emperature at EM(5)0, channel (3)
 rsb=zeros(nr8,4);
 rs=zeros(nr8,1);
 rs2=zeros(nr8,1);
 rst20=zeros(nr8,1);
% (A)erodynamic (R)esistance using (T)emperature at EM(5)0, channel (4)
 rab=zeros(ns,4);
 rat54=zeros(ns,1);
 raa=zeros(ns,1);
 raa2=zeros(ns,1);
 rat20=zeros(ns,1);
 
  % rac,rsc, is the changing resistances see page 203
  rac=zeros(nr8,1);
  rsc=zeros(nr8,1);
 drva=zeros(nr8,1);
 
 
  % (max)mum et rate 
   etm8=zeros(ns,1);
   elt8=zeros(ns,1);
 
   va=zeros(nr8,1);
   va1=zeros(nr8,1);
   va2=zeros(nr8,1);
   advm=zeros(nr8,1);
   adv1=zeros(nr8,1);
   adv2=zeros(nr8,1);
 %  (RA) when soil surface (T)emperature is only (20)
 t20k=20+az;
  for i=1:ns
  if (so(i)==1 && wl(i)==0)
 %  interpolation to DT(8)5G time series
   % accumulative evaporation of 60cm column at dt85g time series
    aet6c8(sst8(i):set8(i))=interp1(dtsd(ssts(i):sets(i)),aet6c(ssts(i):sets(i)),dt8d(sst8(i):set8(i)));
 % Transient evaporation of 60cm column at dt85g time series
    et68(sst8(i):set8(i))=interp1(dtsd(ssts(i):sets(i)),etsgf6(ssts(i):sets(i)),dt8d(sst8(i):set8(i)));
    
    rh538(sst8(i):set8(i))=interp1(dt5d(sst5(i):set5(i)),rh53(sst5(i):set5(i)),dt8d(sst8(i):set8(i)),'spline');  
    rh548(sst8(i):set8(i))=interp1(dt5d(sst5(i):set5(i)),rh54(sst5(i):set5(i)),dt8d(sst8(i):set8(i)),'spline');   
    t538(sst8(i):set8(i))=interp1(dt5d(sst5(i):set5(i)),t53(sst5(i):set5(i)),dt8d(sst8(i):set8(i)));
    t548(sst8(i):set8(i))=interp1(dt5d(sst5(i):set5(i)),t54(sst5(i):set5(i)),dt8d(sst8(i):set8(i)));
    
    
    t548(sst8(i))=t548(sst8(i)+1);
    t538(sst8(i))=t538(sst8(i)+1);
    rh548(sst8(i))=rh548(sst8(i)+1);
    rh538(sst8(i))=rh538(sst8(i)+1);
%      % check the interpolation
%      h=figure;i=4;
%     plot(dt5d(sst5(i):set5(i))-dt5d(sst5(i)),rh53(sst5(i):set5(i)),'r','LineWidth',4);hold on
%     plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),rh538(sst8(i):set8(i)));hold on
%     
%     plot(dt5d(sst5(i):set5(i))-dt5d(sst5(i)),rh54(sst5(i):set5(i)),'r','LineWidth',4);hold on
%     plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),rh548(sst8(i):set8(i)));hold on
%     
%     
%          h=figure;i=4;
%     plot(dt5d(sst5(i):set5(i))-dt5d(sst5(i)),t53(sst5(i):set5(i)),'r','LineWidth',4);hold on
%     plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),t538(sst8(i):set8(i)));hold on
%     
%     plot(dt5d(sst5(i):set5(i))-dt5d(sst5(i)),t54(sst5(i):set5(i)),'r','LineWidth',4);hold on
%     plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),t548(sst8(i):set8(i)));hold on
    
    
    
    
    
    
    
    
%       h=figure;i=4;
%      plot(dtsd(ssts(i):sets(i))-dtsd(ssts(i)),aet6c(ssts(i):sets(i)),'r','LineWidth',4);hold on
%      plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),aet6c8(sst8(i):set8(i)));hold on
    
%  calculating aerodyn% (max)mum et rate 
   [etm8(i),etl8(i)]=max(etsp8(sst8(i):set8(i)));%amic resistance
   
   
   % change etl8 into absolute value
   etl8(i)=etl8(i)+sst8(i);
   
   for k=1:3
      
   
       
       
   % remember, here the relative humidty on the soil surface is not based on the measurement at the sensor, instead, it is based on the kelvin equation.
   rab(i,k)=(func.svp(t538(sst8(i))+az)*1./(t538(sst8(i))+az)-func.svp(tma(k)+az)*rh(k)/(tma(k)+az))*mw/R/rholw./etsp8(sst8(i));    
   
   rsb(sst8(i):set8(i),k)=(func.svp(t538(sst8(i))+az)*1./(t538(sst8(i))+az)-func.svp(tma(k)+az).*rh(k)./(tma(k)+az))*mw/R/rholw./etsp8(sst8(i):set8(i))-rab(i,k);
       
       
   end
   
   
   
   
   k=4; % the fourth session where temperature and humidity is the the one 
   
% remember, here the relative humidty on the soil surface is not based on the measurement at the sensor, instead, it is based on the kelvin equation.
   rab(i,k)=(func.svp(t538(sst8(i))+az)*1./(t538(sst8(i))+az)-func.svp(t548(sst8(i))+az)*rh548(sst8(i))/(t548(sst8(i))+az))*mw/R/rholw./etsp8(sst8(i));    
   rsb(sst8(i):set8(i),k)=(func.svp(t538(sst8(i):set8(i) )+az)*1./(t538(sst8(i):set8(i))+az)-func.svp(t548(sst8(i):set8(i))+az).*rh548(sst8(i))./(t548(sst8(i):set8(i))+az))*mw/R/rholw./etsp8(sst8(i):set8(i))-rab(i,k);
   
   
% %  ra as the first 
%    raa(i)=(func.svp(at54(i)+az)*1/(at54(i)+az)-func.svp(tma+az)*rh/(tma+az))*mw/R/rholw./etsp8(zp8at(i));
% %  et8(zp8at(i)+50) refers the moisture content in the middle (like 0.5)
%    raa2(i)=(func.svp(at54(i)+az)*1/(at54(i)+az)-func.svp(tma+az)*rh/(tma+az))*mw/R/rholw./etsp8(zp8at(i)+50); %et8(zp8at(i)+50) means
% 
%    rat20(i)=(func.svp(t20k)*1/t20k-func.svp(tma+az)*rh/(tma+az))*mw/R/rholw./etsp8(zp8at(i)+50); %et8(zp8at(i)+50) means
%    
%    rat54(i)=(func.svp(t548(zp8at(i)+50)+az)*1/(t548(zp8at(i)+50)+az)-func.svp(tma+az)*rh/(tma+az)) ...
%        *mw/R/rholw./etsp8(zp8at(i)+50); %et8(zp8at(i)+50) means
% %  calculate surface resistance   
%    rs(zp8at(i):set8(i))=(func.svp(at54(i)+az)*1/(at54(i)+az)-func.svp(tma+az)*rh/(tma+az))*mw/R/rholw./etsp8(zp8at(i):set8(i))-raa(i);
% 
%    rs2(zp8at(i):set8(i))=(func.svp(at54(i)+az)*1/(at54(i)+az)-func.svp(tma+az)*rh/(tma+az))*mw/R/rholw./etsp8(zp8at(i):set8(i))-raa2(i);
%    
%    
%    rst20(zp8at(i):set8(i))=(func.svp(t20k)*1/t20k-func.svp(tma+az)*rh/(tma+az))*mw/R/rholw./etsp8(zp8at(i):set8(i))-rat20(i);
%    
%    
%    rst51(zp8at(i):set8(i))=(func.svp(t548(zp8at(i):set8(i))+az)*1./(t548(zp8at(i):set8(i))+az)-func.svp(tma+az)*rh/(tma+az)) ...
%        *mw/R/rholw./etsp8(zp8at(i):set8(i))-rat54(i);
 % calculate the (E)vapora(T)ion from (6)0 column rate based on the relative humidity sensor at the tunnel above the soil surface
 % using DT(8)5G time series (m/s)
 
    et6rh8(sst8(i):set8(i))=func.dv(  (t548(sst8(i):set8(i))+t538(sst8(i):set8(i))+2*az)/2     ).*...
      (func.rhovs(t548(sst8(i):set8(i))+az).*rh548(sst8(i):set8(i))-func.rhovs(t538(sst8(i):set8(i))+az).*rh538(sst8(i):set8(i)))...
      /0.1/rholw;
   
  
  va(sst8(i):set8(i))=-kg/miua/por*( func.svp(t548(sst8(i):set8(i))+az).*rh548(sst8(i):set8(i))- func.svp(60+az)*0.01 )/dz;
  
  va1(sst8(i):set8(i))=-kg/miua/por*( func.svp(t538(sst8(i):set8(i))+az+5).*(rh538(sst8(i):set8(i)))- func.svp(60+az)*0.01 )/dz;
  
  va2(sst8(i):set8(i))=-kg/miua/por*( func.svp(t548(sst8(i):set8(i))+az).*rh548(sst8(i):set8(i))- func.svp(60+az)*0.01 )/dz1;
  
  
  

%    adv1(sst8(i):set8(i))=-va(sst8(i):set8(i)).*func.rhovs(t548(sst8(i):set8(i))+az).*rh548(sst8(i):set8(i))  /rholw*ms2mmd;
%    adv2(sst8(i):set8(i))=-va(sst8(i):set8(i)).*func.rhovs(t538(sst8(i):set8(i))+az).*rh538(sst8(i):set8(i))  /rholw*ms2mmd;
  
  advm(sst8(i):set8(i))=-va(sst8(i):set8(i)).*func.rhovs( (t548(sst8(i):set8(i))+t538(sst8(i):set8(i))+2*az)/2 ).*(rh548(sst8(i):set8(i)) +rh538(sst8(i):set8(i)) )/2/rholw;
  
  adv1(sst8(i):set8(i))=-va1(sst8(i):set8(i)).*func.rhovs( t538(sst8(i):set8(i))+az +5).*(rh538(sst8(i):set8(i)) )/rholw;
  
  adv2(sst8(i):set8(i))=-va2(sst8(i):set8(i)).*func.rhovs( t548(sst8(i):set8(i))+az ).*(rh548(sst8(i):set8(i)) )/rholw;
  
  
  % the difference of water vapor between two sensors used for geting rac
drva(sst8(i):set8(i))=func.rhovs(t538(sst8(i):set8(i))+az).*rh538(sst8(i):set8(i))-func.rhovs(60+az)*0.01;
rac(sst8(i):set8(i))=drva(sst8(i):set8(i))/rholw./adv1(sst8(i):set8(i));
rsc(sst8(i):set8(i))=(func.svp(    t538(sst8(i):set8(i))+az     )*1./(  t538(sst8(i):set8(i))+az   )-func.svp(60+az).*0.01./(60+az)) ...
      *mw/R/rholw./etsp8(sst8(i):set8(i))-rac(sst8(i):set8(i));  
  
  
%    figure;
%     i=4;
%     plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),rsc(sst8(i):set8(i)));hold on
%     
%     
%     plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),rac(sst8(i):set8(i)));hold on
  
  
  
  
  
  end
  
  
  end
        fprintf(1, strcat('Interpolation finished \n'));

%% plot mean temperature against depth at each stable conditions

    figure
    semilogy(mc8f(sst8(5):iv:set8(5)),rsb(sst8(5):iv:set8(5),4),'ro');hold on    %M14 96.19
    semilogy(mc8f(sst8(4):iv:set8(4)),rsb(sst8(4):iv:set8(4),4),'b.');  
%     xlim([0 1])
%    h=figure;
%    set(gcf,'Units','normalized', 'WindowStyle','docked','OuterPosition',[0 0 1 1]); %maximized the window
%    subplot('Position',[0.05 0.65 0.18 0.32])
%    
% %  h1=plot([at53(1),at55(1),at54(1),at83(1),at84(1)],[es52,es55,es51,es83,es84],'r+','LineWidth',lw);hold on
% %    i=2;h2=plot([at53(i),at54(i),at82(i),at81(i)],[es53,es54,es82,es81],'go','LineWidth',lw,'MarkerSize',ms);hold on
% %    i=3;h3=plot([at53(i),at54(i),at82(i),at81(i)],[es53,es54,es82,es81],'b*','LineWidth',lw,'MarkerSize',ms);hold on
%    i=4;h4=plot([at53(i),at54(i),at82(i),at81(i)],[es53,es54,es82,es81],'cx','LineWidth',lw,'MarkerSize',ms);hold on
%    i=5;h5=plot([at53(i),at54(i),at82(i),at81(i)],[es53,es54,es82,es81],'ms','LineWidth',lw,'MarkerSize',ms);hold on
% %    i=6;h6=plot([at53(i),at54(i),at82(i),at81(i)],[es53,es54,es82,es81],'yd','LineWidth',lw,'MarkerSize',ms);hold on
% %    i=7;h7=plot([at53(i),at54(i),at82(i),at81(i)],[es53,es54,es82,es81],'k^','LineWidth',lw,'MarkerSize',ms);hold on
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
% %   line([atb(1) atb(1)],yl,'LineStyle',':','Color','r','LineWidth',lw);
% %   line([atb(2) atb(2)],yl,'LineStyle',':','Color','g','LineWidth',lw);
% %   line([atb(3) atb(3)],yl,'LineStyle',':','Color','b','LineWidth',lw);
% %   line([atb(4) atb(4)],yl,'LineStyle',':','Color','c','LineWidth',lw);
% %   line([atb(5) atb(5)],yl,'LineStyle',':','Color','m','LineWidth',lw);
% %   line([atb(6) atb(6)],yl,'LineStyle',':','Color','y','LineWidth',lw);
% %   line([atb(7) atb(7)],yl,'LineStyle',':','Color','k','LineWidth',lw);
% %   line([atb(8) atb(8)],yl,'LineStyle','-.','Color','r','LineWidth',lw);
% %   line([atb(9) atb(9)],yl,'LineStyle','-.','Color','g','LineWidth',lw);
% %   line([atb(10) atb(10)],yl,'LineStyle','-.','Color','b','LineWidth',lw);
%     
%  ax1=gca;set(ax1,'FontSize',tfs)  
%  xlabel('Temperature ( ^\circ C)','FontSize',fs,'FontWeight','bold')
%  ylabel('Depth (m)','FontSize',fs,'FontWeight','bold')
% %   hleg1 = legend([h4,h5],[cs(4),cs(5)],'Location','SouthEast');
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
% % %   line([arhb(1) arhb(1)],yl,'LineStyle',':','Color','r','LineWidth',lw);
% % %   line([arhb(2) arhb(2)],yl,'LineStyle',':','Color','g','LineWidth',lw);
% % %   line([arhb(3) arhb(3)],yl,'LineStyle',':','Color','b','LineWidth',lw);
% % %   line([arhb(4) arhb(4)],yl,'LineStyle',':','Color','c','LineWidth',lw);
% % %   line([arhb(5) arhb(5)],yl,'LineStyle',':','Color','m','LineWidth',lw);
% % %   line([arhb(6) arhb(6)],yl,'LineStyle',':','Color','y','LineWidth',lw);
% % %   line([arhb(7) arhb(7)],yl,'LineStyle',':','Color','k','LineWidth',lw);
% % %   line([arhb(8) arhb(8)],yl,'LineStyle','-.','Color','r','LineWidth',lw);
% % %   line([arhb(9) arhb(9)],yl,'LineStyle','-.','Color','g','LineWidth',lw);
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
%  saveas(h,['Medium_60Analysis_Transient_60_Temperature','fig'],'fig')
%   fprintf(1, strcat('First Graph finished \n'));
%  h=figure;
%  set(gcf,'Units','normalized', 'WindowStyle','docked','OuterPosition',[0 0 1 1]); %maximized the window
% %% plot Moisture content after translation
% subplot('Position',[0.55 0.70 0.43 0.28])
% % i=3;h3=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),mc8f(zp8at(i):set8(i)),'g','LineWidth',lw);hold on
% i=4;h4=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),mc8f(zp8at(i):set8(i)),'c','LineWidth',lw);hold on
% i=5;h5=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),mc8f(zp8at(i):set8(i)),'m','LineWidth',lw);hold on
% % i=6;h6=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),mc8f(zp8at(i):set8(i)),'y','LineWidth',lw);hold on
% % i=7;h7=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),mc8f(zp8at(i):set8(i)),'k','LineWidth',lw);hold on
% axis([0,8,0,1])
% xlabel('Time (day)','FontSize',fs)
% ylabel('Sautration (-)','FontSize',fs)
% %% plot accumulative evaporation rate after translation
% subplot('Position',[0.55 0.38 0.43 0.28])
% % i=3;h3=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),aet6c(zpsat(i):sets(i))-aet6c(zpsat(i)),'g','LineWidth',lw);hold on
% i=4;h4=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),aet6c(zpsat(i):sets(i))-aet6c(zpsat(i)),'c','LineWidth',lw);hold on
% i=5;h5=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),aet6c(zpsat(i):sets(i))-aet6c(zpsat(i)),'m','LineWidth',lw);hold on
% % i=6;h6=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),aet6c(zpsat(i):sets(i))-aet6c(zpsat(i)),'y','LineWidth',lw);hold on
% % i=7;h7=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),aet6c(zpsat(i):sets(i))-aet6c(zpsat(i)),'k','LineWidth',lw);hold on
% 
% ylabel('AET (mm)','FontSize',fs)
% xlabel('Time (day)','FontSize',fs)
% %% Plot Transient evaportion rate after translation as a function of time
% subplot('Position',[0.55 0.05 0.43 0.28])
% % i=3;h2=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf6(zpsat(i):sets(i))*ms2mmd,'g','LineWidth',lw);hold on
% % i=4;h4=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf6(zpsat(i):sets(i))*ms2mmd,'c','LineWidth',lw);hold on
% % i=5;h5=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf6(zpsat(i):sets(i))*ms2mmd,'m','LineWidth',lw);hold on
% % i=6;h6=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf6(zpsat(i):sets(i))*ms2mmd,'y','LineWidth',lw);hold on
% % i=7;h7=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf6(zpsat(i):sets(i))*ms2mmd,'k','LineWidth',lw);hold on
% 
% for i=1:ns
%  if (wl(i)==0 && so(i)==1)
% % plot(dtsd(ssts(i):sets(i))-dtsd(ssts(i)),etsgf9(ssts(i):sets(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),etsp8(sst8(i):set8(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% %plot(dt8d(
% 8(i):iv:set8(i))-dt8d(sst8(i)),advm(sst8(i):iv:set8(i))*ms2mmd,'+','LineWidth',lw,'color',col(i,:));hold on
% plot(dt8d(sst8(i):iv:set8(i))-dt8d(sst8(i)),adv1(sst8(i):iv:set8(i))*ms2mmd,'o','LineWidth',lw,'color',col(i,:));hold on
% %plot(dt8d(sst8(i):iv:set8(i))-dt8d(sst8(i)),adv2(sst8(i):iv:set8(i))*ms2mmd,'*','LineWidth',lw,'color',col(i,:));hold on
% 
%  end
% end
% 
% 
% xlabel('Time (day)','FontSize',fs)
% ylabel('ET rate (mm/day)','FontSize',fs)
% %% Plot Transient evaportion rate as a function of accumulative evaporation rate
% subplot('Position',[0.04 0.68 0.43 0.28])
% % i=3;h2=plot(aet6c(zpsat(i):sets(i))-aet6c(zpsat(i)),etsgf6(zpsat(i):sets(i))*ms2mmd,'g','LineWidth',lw);hold on
% i=4;h4=plot(aet6c(zpsat(i):sets(i))-aet6c(zpsat(i)),etsgf6(zpsat(i):sets(i))*ms2mmd,'c','LineWidth',lw);hold on
% i=5;h5=plot(aet6c(zpsat(i):sets(i))-aet6c(zpsat(i)),etsgf6(zpsat(i):sets(i))*ms2mmd,'m','LineWidth',lw);hold on
% % i=6;h6=plot(aet6c(zpsat(i):sets(i))-aet6c(zpsat(i)),etsgf6(zpsat(i):sets(i))*ms2mmd,'y','LineWidth',lw);hold on
% % i=7;h7=plot(aet6c(zpsat(i):sets(i))-aet6c(zpsat(i)),etsgf6(zpsat(i):sets(i))*ms2mmd,'k','LineWidth',lw);hold on
% 
% xlabel('AET (mm)','FontSize',fs)
% ylabel('ET rate (mm/day)','FontSize',fs)
% %% Plot Transient evaporation rate as a function of saturation
% subplot('Position',[0.04 0.36 0.43 0.28])
% % i=3;h2=plot(mc8f(zp8at(i):set8(i)),etsp8(zp8at(i):set8(i))*ms2mmd,'g','LineWidth',lw);hold on
% i=4;h4=plot(mc8f(zp8at(i):set8(i)),etsp8(zp8at(i):set8(i))*ms2mmd,'c','LineWidth',lw);hold on
% i=5;h5=plot(mc8f(zp8at(i):set8(i)),etsp8(zp8at(i):set8(i))*ms2mmd,'m','LineWidth',lw);hold on
% % i=6;h6=plot(mc8f(zp8at(i):set8(i)),etsp8(zp8at(i):set8(i))*ms2mmd,'y','LineWidth',lw);hold on
% % i=7;h7=plot(mc8f(zp8at(i):set8(i)),etsp8(zp8at(i):set8(i))*ms2mmd,'k','LineWidth',lw);hold on
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
% 
% 
% % i=3;h2=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),et6rh8(zp8at(i):set8(i))*ms2mmd,'g','LineWidth',lw);hold on
% i=4;h4=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),et6rh8(zp8at(i):set8(i))*ms2mmd,'c','LineWidth',lw);hold on
% i=5;h5=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),et6rh8(zp8at(i):set8(i))*ms2mmd,'m','LineWidth',lw);hold on
% % i=6;h6=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),et6rh8(zp8at(i):set8(i))*ms2mmd,'y','LineWidth',lw);hold on
% % i=7;h7=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),et6rh8(zp8at(i):set8(i))*ms2mmd,'k','LineWidth',lw);hold on
% 
% xlabel('Time (day)','FontSize',fs)
% ylabel('ET by RH (mm/day)','FontSize',fs)
% axis([0,8,0,1])
% 
% 
% saveas(h,['Medium_60_Analysis_Transient_60_ETrate','.fig'],'fig')
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
% 
% % i=3;h3=plot(aet6c(zpsat(i):sets(i))-aet6c(zpsat(i)),etsgf6(zpsat(i):sets(i))/etsgf6(zpsat(i)),'g','LineWidth',lw);hold on
% i=4;h4=plot(aet6c(zpsat(i):sets(i))-aet6c(zpsat(i)),etsgf6(zpsat(i):sets(i))/etsgf6(zpsat(i)),'c','LineWidth',lw);hold on
% i=5;h5=plot(aet6c(zpsat(i):sets(i))-aet6c(zpsat(i)),etsgf6(zpsat(i):sets(i))/etsgf6(zpsat(i)),'m','LineWidth',lw);hold on
% % i=6;h6=plot(aet6c(zpsat(i):sets(i))-aet6c(zpsat(i)),etsgf6(zpsat(i):sets(i))/etsgf6(zpsat(i)),'y','LineWidth',lw);hold on
% % i=7;h7=plot(aet6c(zpsat(i):sets(i))-aet6c(zpsat(i)),etsgf6(zpsat(i):sets(i))/etsgf6(zpsat(i)),'k','LineWidth',lw);hold on
% 
% xlabel('AET (mm)','FontSize',fs)
% ylabel('ET rate (-)','FontSize',fs)
% 
% 
% %% Plot non-dimensional evaporation rate as a function of saturation
% 
% subplot('Position',[0.04 0.36 0.43 0.28])
% % i=3;h3=plot(mc8f(zp8at(i):set8(i)),etsp8(zp8at(i):set8(i))/etsp8(zp8at(i)),'g','LineWidth',lw);hold on
% i=4;h4=plot(mc8f(zp8at(i):set8(i)),etsp8(zp8at(i):set8(i))/etsp8(zp8at(i)),'c','LineWidth',lw);hold on
% i=5;h5=plot(mc8f(zp8at(i):set8(i)),etsp8(zp8at(i):set8(i))/etsp8(zp8at(i)),'m','LineWidth',lw);hold on
% % i=6;h6=plot(mc8f(zp8at(i):set8(i)),etsp8(zp8at(i):set8(i))/etsp8(zp8at(i)),'y','LineWidth',lw);hold on
% % i=7;h7=plot(mc8f(zp8at(i):set8(i)),etsp8(zp8at(i):set8(i))/etsp8(zp8at(i)),'k','LineWidth',lw);hold on
% 
% xlabel('Saturation (-)','FontSize',fs)
% ylabel('ET rate (-)','FontSize',fs)
% axis([0,1,0,1.2])
% %% Plot non-dimensional evaporation rate calculated by relative humidity against time
% subplot('Position',[0.04 0.04 0.43 0.28])
% 
% % i=3;h3=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),et6rh8(zp8at(i):set8(i))/et6rh8(zp8at(i)),'g','LineWidth',lw);hold on
% i=4;h4=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),et6rh8(zp8at(i):set8(i))/et6rh8(zp8at(i)),'c','LineWidth',lw);hold on
% i=5;h5=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),et6rh8(zp8at(i):set8(i))/et6rh8(zp8at(i)),'m','LineWidth',lw);hold on
% % i=6;h6=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),et6rh8(zp8at(i):set8(i))/et6rh8(zp8at(i)),'y','LineWidth',lw);hold on
% % i=7;h6=plot(dt8d(zp8at(i):set8(i))-dt8d(zp8at(i)),et6rh8(zp8at(i):set8(i))/et6rh8(zp8at(i)),'k','LineWidth',lw);hold on
% 
% xlabel('Time (day)','FontSize',fs)
% ylabel('ET by RH (-)','FontSize',fs)
% axis([0,8,0,1])
% 
% 
% %% Plot non-dimensional evaporation rate calculated by relative humidity against time
% subplot('Position',[0.55 0.36 0.43 0.28])
% 
% % i=3;h3=plot(aet6c8(zp8at(i):set8(i))-aet6c8(zp8at(i)),et6rh8(zp8at(i):set8(i))/et6rh8(zp8at(i)),'g','LineWidth',lw);hold on
% i=4;h4=plot(aet6c8(zp8at(i):set8(i))-aet6c8(zp8at(i)),et6rh8(zp8at(i):set8(i))/et6rh8(zp8at(i)),'c','LineWidth',lw);hold on
% i=5;h5=plot(aet6c8(zp8at(i):set8(i))-aet6c8(zp8at(i)),et6rh8(zp8at(i):set8(i))/et6rh8(zp8at(i)),'m','LineWidth',lw);hold on
% % i=6;h6=plot(aet6c8(zp8at(i):set8(i))-aet6c8(zp8at(i)),et6rh8(zp8at(i):set8(i))/et6rh8(zp8at(i)),'y','LineWidth',lw);hold on
% % i=7;h6=plot(aet6c8(zp8at(i):set8(i))-aet6c8(zp8at(i)),et6rh8(zp8at(i):set8(i))/et6rh8(zp8at(i)),'k','LineWidth',lw);hold on
% 
% xlabel('Time (day)','FontSize',fs)
% ylabel('ET by RH (-)','FontSize',fs)
% 
% 
% %% Plot Transient evaportion rate after translation as a function of time
% subplot('Position',[0.55 0.05 0.43 0.28])
% % i=3;h3=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf6(zpsat(i):sets(i))/etsgf6(zpsat(i)),'g','LineWidth',lw);hold on
% %i=4;h4=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf6(zpsat(i):sets(i))/etsgf6(zpsat(i)),'c','LineWidth',lw);hold on
% %i=5;h5=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf6(zpsat(i):sets(i))/etsgf6(zpsat(i)),'m','LineWidth',lw);hold on
% % i=6;h6=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf6(zpsat(i):sets(i))/etsgf6(zpsat(i)),'y','LineWidth',lw);hold on
% % i=7;h7=plot(dtsd(zpsat(i):sets(i))-dtsd(zpsat(i)),etsgf6(zpsat(i):sets(i))/etsgf6(zpsat(i)),'k','LineWidth',lw);hold on
% 
% 
% for i=1:ns
%  if (wl(i)==0 && so(i)==1)
% % plot(dtsd(ssts(i):sets(i))-dtsd(ssts(i)),etsgf9(ssts(i):sets(i))*ms2mmd,'-','LineWidth',lw,'color',col(i,:));hold on
% plot(dt8d(sst8(i):set8(i))-dt8d(sst8(i)),etsp8(sst8(i):set8(i))/etsp8(sst8(i)),'-','LineWidth',lw,'color',col(i,:));hold on
% %plot(dt8d(sst8(i):iv:set8(i))-dt8d(sst8(i)),advm(sst8(i):iv:set8(i))/advm(sst8(i)),'+','LineWidth',lw,'color',col(i,:));hold on
% plot(dt8d(sst8(i):iv:set8(i))-dt8d(sst8(i)),adv1(sst8(i):iv:set8(i))/adv1(sst8(i)),'o','LineWidth',lw,'color',col(i,:));hold on
% %plot(dt8d(sst8(i):iv:set8(i))-dt8d(sst8(i)),adv2(sst8(i):iv:set8(i))/adv2(sst8(i)),'*','LineWidth',lw,'color',col(i,:));hold on
% 
%  end
% end
% 
% 
% xlabel('AET (mm)','FontSize',fs)
% ylabel('non-ET rate (-)','FontSize',fs)
% 
% 
% saveas(h,['Medium_60_Analysis_Transient_60_nondimensional','.fig'],'fig')
% fprintf(1, strcat('non-dimensional Graph finished \n'));
% % %% Plot Surface resistance as a function of saturation
% % h=figure;
% % set(gcf,'Units','normalized', 'WindowStyle','docked','OuterPosition',[0 0 1 1]); %maximized the window
% % 
% %  subplot('Position',[0.04 0.70 0.43 0.28])
% %   ha=plot(0.01:0.05:1,func.rs1994(0.01:0.05:1,0.375,por),'k-','LineWidth',lw);hold on % van der grivend 1994
% %   hb=plot(0.03:0.01:0.45,func.rs1986(0.03:0.01:0.45,por),'k-.','LineWidth',lw);hold on % camilo 1986
% %   hc=plot(0.01:0.05:1,func.rs1996(0.01:0.05:1,por),':','LineWidth',lw);hold on % Daamen 1996
% %   hd=plot(0.01:0.05:1,func.rs1984(0.01:0.05:1),'--','LineWidth',lw);hold on % Sun 1984
% %   
% %   
% % %   i=3;h3=plot(mc8f(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'g','LineWidth',lw);hold on
% %   i=4;h4=plot(mc8f(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'c','LineWidth',lw);hold on
% %   i=5;h5=plot(mc8f(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'m','LineWidth',lw);hold on
% % %   i=6;h6=plot(mc8f(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'y','LineWidth',lw);hold on
% % %   i=7;h7=plot(mc8f(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'k','LineWidth',lw);hold on
% %   
% %   
% % 
% %   
% %   axis([0,1,0,600])
% %   hleg1 = legend([ha,hb,hc,hd],{sprintf('Van der Griend et al. [1994]'),sprintf('Camillo et al. [1986]'),sprintf('Daamen et al. [1996]'), ...
% %   sprintf('Sun [1984]')},'Location','NorthEast');
% %   xlabel('Saturation (-)','FontSize',fs)
% %   ylabel('rs (s/m)','FontSize',fs)
% % 
% % %% Plot surface resistance on log profile
% % subplot('Position',[0.04 0.36 0.43 0.28])
% %   ha=semilogy(0.01:0.05:1,func.rs1994(0.01:0.05:1,0.375,por),'k-','LineWidth',lw);hold on % van der grivend 1994
% %   hb=semilogy(0.03:0.01:0.45,func.rs1986(0.03:0.01:0.45,por),'k-.','LineWidth',lw);hold on % camilo 1986
% %   hc=semilogy(0.01:0.05:1,func.rs1996(0.01:0.05:1,por),':','LineWidth',lw);hold on % Daamen 1996
% %   hd=semilogy(0.01:0.05:1,func.rs1984(0.01:0.05:1),'--','LineWidth',lw);hold on % Sun 1984
% %   
% %   
% %   
% % %   i=3;h3=semilogy(mc8f(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'g','LineWidth',lw);hold on
% %   i=4;h4=semilogy(mc8f(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'c','LineWidth',lw);hold on
% %   i=5;h5=semilogy(mc8f(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'m','LineWidth',lw);hold on
% % %   i=6;h6=semilogy(mc8f(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'y','LineWidth',lw);hold on
% % %   i=7;h7=semilogy(mc8f(zp8at(i):set8(i)),rs(zp8at(i):set8(i)),'k','LineWidth',lw);hold on
% %   
% %   axis([0,1,1,1000])
% %   hleg1 = legend([ha,hb,hc,hd],{sprintf('Van der Griend et al. [1994]'),sprintf('Camillo et al. [1986]'),sprintf('Daamen et al. [1996]'), ...
% %   sprintf('Sun [1984]')},'Location','NorthEast');
% %   xlabel('Saturation (-)','FontSize',fs)
% %   ylabel('rs (s/m)','FontSize',fs)
% % %% Plot Surface resistance as a function of saturation with translation
% %  subplot('Position',[0.04 0.04 0.43 0.3])
% %   ha=plot(0.01:0.05:1,func.rs1994(0.01:0.05:1,0.375,por),'k-','LineWidth',lw);hold on % van der grivend 1994
% %   hb=plot(0.03:0.01:0.45,func.rs1986(0.03:0.01:0.45,por),'k-.','LineWidth',lw);hold on % camilo 1986
% %   hc=plot(0.01:0.05:1,func.rs1996(0.01:0.05:1,por),':','LineWidth',lw);hold on % Daamen 1996
% %   hd=plot(0.01:0.05:1,func.rs1984(0.01:0.05:1),'--','LineWidth',lw);hold on % Sun 1984
% % %   i=3;h3=plot(mc8f(zp8at(i):set8(i)),rs2(zp8at(i):set8(i)),'g','LineWidth',lw);hold on
% %   i=4;h4=plot(mc8f(zp8at(i):set8(i)),rs2(zp8at(i):set8(i)),'c','LineWidth',lw);hold on
% %   i=5;h5=plot(mc8f(zp8at(i):set8(i)),rs2(zp8at(i):set8(i)),'m','LineWidth',lw);hold on
% % %   i=6;h6=plot(mc8f(zp8at(i):set8(i)),rs2(zp8at(i):set8(i)),'y','LineWidth',lw);hold on
% % %   i=7;h7=plot(mc8f(zp8at(i):set8(i)),rs2(zp8at(i):set8(i)),'k','LineWidth',lw);hold on
% %   
% %   
% % %   i=3;g3=plot(mc8f(zp8at(i):set8(i)),rst20(zp8at(i):set8(i)),'g-.','LineWidth',lw);hold on
% %   i=4;g4=plot(mc8f(zp8at(i):set8(i)),rst20(zp8at(i):set8(i)),'c-.','LineWidth',lw);hold on
% %   i=5;g5=plot(mc8f(zp8at(i):set8(i)),rst20(zp8at(i):set8(i)),'m-.','LineWidth',lw);hold on
% % %   i=6;g6=plot(mc8f(zp8at(i):set8(i)),rst20(zp8at(i):set8(i)),'y-.','LineWidth',lw);hold on
% % %   i=7;g7=plot(mc8f(zp8at(i):set8(i)),rst20(zp8at(i):set8(i)),'k-.','LineWidth',lw);hold on
% %   
% %   
% % %   i=3;g32=plot(mc8f(zp8at(i):set8(i)),rst51(zp8at(i):set8(i)),'g--','LineWidth',lw);hold on
% % %   i=4;g42=plot(mc8f(zp8at(i):set8(i)),rst51(zp8at(i):set8(i)),'c--','LineWidth',lw);hold on
% % %   i=5;g52=plot(mc8f(zp8at(i):set8(i)),rst51(zp8at(i):set8(i)),'m--','LineWidth',lw);hold on
% % %   i=6;g62=plot(mc8f(zp8at(i):set8(i)),rst51(zp8at(i):set8(i)),'y--','LineWidth',lw);hold on
% % %   i=7;g72=plot(mc8f(zp8at(i):set8(i)),rst51(zp8at(i):set8(i)),'k--','LineWidth',lw);hold on
% %   
% %   
% %   
% %   axis([0,1,-100,600])
% %   hleg1 = legend([ha,hb,hc,hd],{sprintf('Van der Griend et al. [1994]'),sprintf('Camillo et al. [1986]'),sprintf('Daamen et al. [1996]'), ...
% %     sprintf('Sun [1984]')},'Location','NorthEast');
% %   xlabel('Saturation (-)','FontSize',fs)
% %   ylabel('rs (s/m)','FontSize',fs)
% %  %% save file
% %  saveas(h,['Medium_60_Analysis_Transient_60_Resistance','.fig'],'fig')
%  fprintf(1, strcat('Analysis finished \n'));
